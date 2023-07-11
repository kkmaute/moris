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

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "GEN_typedefs.hpp"
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
                                      << " Ii: " << Ii << std::endl;
                            tInterface_21 = tInterface_21 + "dbl_iside_p0_" + std::to_string( Ik ) + "_p1_" + std::to_string( Ii ) + ",";
                            std::cout << tInterface_21 << Ii << std::endl;
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
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // Inlet velocity function
    void
    Func_Inlet_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
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
    get_phase_index( const moris::ge::Geometry_Bitset& aGeometrySigns )
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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "39,39,39" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "4.0,4.0,4.0" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-2.0,-2.0,-2.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4,5,6" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    void
    XTKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "multigrid", isMultigrid );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
    }

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "intersection_mode", "LEVEL_SET" );
        tParameterlist( 0 )( 0 ).set( "number_of_phases", (sint)( tNumBulkPhases + 1 ) );
        tParameterlist( 0 )( 0 ).set( "phase_function_name", "get_phase_index" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", "Voxel_Grain_Structure_Gen.exo" );

        // Geometry parameter lists
        tParameterlist( 1 )( 0 ) = prm::create_voxel_field_parameter_list();

        std::string tMORISROOT = moris::get_base_moris_dir();
        std::string tPath      = tMORISROOT + "projects/EXA/thermal/diffusion/Voxel_Grain_Structure/Voxel_file";

        tParameterlist( 1 )( 0 ).set( "voxel_field_file", tPath );
        tParameterlist( 1 )( 0 ).set( "domain_dimensions", "4.0,4.0,4.0" );
        tParameterlist( 1 )( 0 ).set( "domain_offset", "-2.0,-2.0,-2.0" );
        tParameterlist( 1 )( 0 ).set( "number_of_refinements", "0" );
        tParameterlist( 1 )( 0 ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( 0 ).set( "discretization_mesh_index", -1 );
    }

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        Name_func();
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity_1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 3
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity_2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInletTemp" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 6
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSideFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_2" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_1,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_3" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model 2
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusion_4" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity_2,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPDirichletNitscheT_1" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity_1,Material" );
        tSPCounter++;

        // create parameter list for stabilization parameter 2
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPDirichletNitscheT_2" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity_2,Material" );
        tSPCounter++;

        // create parameter list for stabilization parameter 8
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp_1" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.005" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity_1,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp_2" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.005" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity_2,Material" );
        tSPCounter++;

        // Temperature - Shell - PCM
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitsche_11" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity_1,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitsche_12" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity_1,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity_2,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitsche_22" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity_2,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPInterfaceNitsche_21" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropConductivity_2,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropConductivity_1,Material" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulk_1 );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionBulk_2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulk_2 );
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInletTemp_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_1,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPDirichletNitscheT_1,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tRightBC_1 );
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGInletTemp_2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropInletTemp,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusion_3,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPDirichletNitscheT_2,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tRightBC_2 );
        tIWGCounter++;

        // create parameter list for IWG 11
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGCylinderFluxTemp" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropSideFlux,Neumann" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLeftBC );
        tIWGCounter++;

        // Temperature - Shell - Shell
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_1" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMDiffusion_2,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_11     ,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface_11 );
        tIWGCounter++;

        // Temperature - Shell - Shell
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_2" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_12     ,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface_12 );
        tIWGCounter++;

        // Temperature - Shell - Shell
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_3" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMDiffusion_4,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_22     ,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface_22 );
        tIWGCounter++;

        // Temperature - Shell - Shell
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGInterfaceShellShellTEMP_4" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", std::string( "TEMP" ) );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMDiffusion_3,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMDiffusion_1,Diffusion" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPInterfaceNitsche_21     ,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterface_21 );
        tIWGCounter++;

        if ( isGhost )
        {
            // create parameter list for IWG 16
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "GPenTemp_1" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp_1,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tGhost_1 );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "GPenTemp_2" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp_2,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tGhost_2 );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 3
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulk );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 7 );
        for ( uint Ik = 0; Ik < 7; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-04 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 1 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "order_adofs_by_host", false );
        tParameterlist( 0 )( 0 ).set( "multigrid", isMultigrid );
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Voxel_Grain_Structure.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tBulk );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "TEMP" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTEMP" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
