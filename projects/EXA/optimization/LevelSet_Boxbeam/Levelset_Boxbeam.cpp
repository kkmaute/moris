/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Levelset_Boxbeam.cpp
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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    //------------------------------------------------------------------------------
    // Main problem parameters

    real tInitialStrainEnergy = 3.26439;
    real tInitialPerimeter    = 29.0669;
    real tPerimeterPenalty    = 0.2;

    real tMaxMass = 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.01;
    int  tMMAMaxIter  = 2;

    int         tInitialRef    = 2;
    std::string tLagrangeOrder = "1";
    std::string tBsplineOrder  = "1,1";

    real tElementEdgeLength = 1.0 / 15.0 / pow( 2, tInitialRef );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    int  tNumHoleX     = 15;
    int  tNumHoleY     = 5;
    real tHoleDiameter = 0.065;

    real tWallThickness = 0.05;

    real tBsplineLimit = 0.06;

    bool tIsOpt    = true;
    bool tUseGhost = true;

    std::string tName = "Levelset_Boxbeam";

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = tName + ".so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets    = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost    = "ghost_p2";
    std::string tInteriorGhost = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets =
            "SideSet_1_n_p2,SideSet_1_c_p2"
            "SideSet_3_n_p2,SideSet_3_c_p2";

    std::string tInterfaceVoidSSets = "iside_b0_2_b1_0,iside_b0_1_b1_0";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    //------------------------------------------------------------------------------
    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj3 = tPerimeterPenalty * aCriteria( 3 ) / tInitialPerimeter;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        std::cout << "min ADV                  = " << aADVs.min() << " \n";
        std::cout << "max ADV                  = " << aADVs.max() << " \n"
                  << std::flush;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 2 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 2 ) << " \n";
        std::cout << "Constraint = " << tConstraints( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenalty / tInitialPerimeter;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 2 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < tLoadLimitY )
        {
            aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
        }
        else
        {
            aPropMatrix = { { 0.0 }, { 0.0 } };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOpt );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
        tParameterlist( 0 )( 0 ).set( "restart_file", "" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_gcmma_parameter_list();
        tParameterlist( 2 )( 0 ).set( "step_size", tMMAStepSize );
        tParameterlist( 2 )( 0 ).set( "penalty", tMMAPenalty );
        tParameterlist( 2 )( 0 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "45,15" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "3,1" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "0.0,0.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", std::to_string( tInitialRef ) );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", false );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {

        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        tParameterlist( 0 )( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        tParameterlist( 0 )( 0 ).set( "print_phase_table", true );

        // init geometry counter
        uint tGeoCounter = 0;

        // outer frame
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 1.5, 0.5, 1.5 - tWallThickness, 0.5 - tWallThickness, 24.0, 0.1, 0.01, 0.0 );
        tGeoCounter++;

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 0.0, 0.0, tHoleDiameter, tHoleDiameter, 8.0, std::pow( tHoleDiameter * tHoleDiameter, 0.5 ), 0.1, 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "lower_bound_x", 0.12 );
        tParameterlist( 1 )( tGeoCounter ).set( "upper_bound_x", 2.88 );
        tParameterlist( 1 )( tGeoCounter ).set( "lower_bound_y", 0.12 );
        tParameterlist( 1 )( tGeoCounter ).set( "upper_bound_y", 0.88 );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_fields_x", tNumHoleX );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_fields_y", tNumHoleY );

        if ( tIsOpt )
        {
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 0 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -tBsplineLimit );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", tBsplineLimit );
            tGeoCounter++;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0e-6" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropTraction" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Neumann_U" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Frame" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Interior" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheFrameInteriorInterface" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhost_Frame" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhost_Interior" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Frame" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Frame" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameSupportSSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameLoadSSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGFrameInteriorInterface" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheFrameInteriorInterface,NitscheInterface" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameInteriorDSets );
        tIWGCounter++;

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGhostFrame" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhost_Frame,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGhostInterior" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhost_Interior,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost );
            tIWGCounter++;
        }
        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInteriorSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Frame" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Interior" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInteriorSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInterfaceVoidSSets );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1.00 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.00 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 1 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", 1 );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", 1.0 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "UX", 1 );
        tParameterlist( 0 )( 0 ).set( "UY", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomainSets );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkUX,IQIBulkUY" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
