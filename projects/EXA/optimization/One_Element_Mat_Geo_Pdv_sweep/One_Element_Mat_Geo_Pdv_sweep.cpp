/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * One_Element_Mat_Geo_Pdv_sweep.cpp
 *
 */

#include <string>
#include <iostream>
#include <sstream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{

    std::string tName = "One_Element_Mat_Geo_Pdv_sweep";

    std::string tLibraryName   = tName + ".so";
    std::string tGENOutputFile = "GEN_" + tName + ".exo";

    std::string tOutputFileName = tName + ".exo";

    bool tIs3D             = false;
    bool tIsOpt            = true;
    bool tUseGhost         = false;
    bool tUseAbsoluteValue = true;

    //-------------------------------
    // Opt constant_parameters

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.021;
    int  tMMAMaxIter  = 400;

    real tInitialStrainEnergy = 1437.01;
    real tMaxMass             = 40;

    real tPerimeterPenalty = 0.05;
    real tInitialPerimeter = 60.0;

    real tRegularizationPenalty = 0.05;    // 0.02
    real tInitialRegularization = 60;

    //-------------------------------

    real phi_sh = 0.5;
    real phi_rt = 3.0;

    int tDispOrder = 1;

    real tBsplineLimitTop    = 1.0;
    real tBsplineLimitBottom = 0.0;

    real tElementEdgeLength = 1.0;

    //------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }
    //-------------------------------
    std::string tBulkSets = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tVoidSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tDirichletSets = "SideSet_4_n_p0,SideSet_4_c_p0";
    std::string tLoadSets      = "SideSet_2_n_p1,SideSet_2_c_p1";

    std::string tInterfaceVoidSets = "iside_b0_1_b1_0";

    std::string tVoidInterfaceSets = "iside_b0_0_b1_1";

    std::string tInteriorGhost = "ghost_p0";

    std::string tTotalDomain = tBulkSets;

    std::string tTotalDomainAGhost = tTotalDomain;

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 20.0 )
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
            }
        }
        else
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { 0.0 }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { 0.0 } };
            }
        }
    }

    real
    Const_Geometry(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        return norm( aCoordinates ) - 3.0;
    }

    //------------------------------------------------------------------------------
    void
    tYoungsFunc( moris::Matrix< moris::DDRMat >&      aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix = tHCT * std::pow( tDensity, tBeta );
    }

    //------------------------------------------------------------------------------
    void
    tDerYoungsFunc( moris::Matrix< moris::DDRMat >&   aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        // if( tDensity <0 ){ tDensity = 0.0001; }
        // if( tDensity >1 ){ tDensity = 1; }

        // FIXME density shift missing

        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() * tBeta * tHCT( 0 ) * std::pow( tDensity, tBeta - 1 ) / ( 1 - phi_sh );
        // aPropMatrix = tBeta * tHCT * std::pow( tDensity, tBeta -1 ) / ( 1 - phi_sh );
    }

    //------------------------------------------------------------------------------

    void
    tDensityFunc( moris::Matrix< moris::DDRMat >&     aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix.set_size( 1, 1, tDensity );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc( moris::Matrix< moris::DDRMat >&  aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() / ( 1 - phi_sh );
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

        tObjectives( 0, 0 ) = obj1;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy            = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";

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
        tConstraints( 0 ) = aCriteria( 1 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 1 ) << " \n";
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

        tDConstraintDCriteria( 1 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );
        aParameterLists.set( "reinitialize_interface_iter", 10000 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        aParameterLists.set( "hdf5_path", "shape_opt_test.hdf5" );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", "1e-6" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", "1,1" );
        aParameterLists.set( "domain_dimensions", "10,10" );
        aParameterLists.set( "domain_offset", "-5.0,-5.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1,1" );
        aParameterLists.set( "lagrange_pattern", "0,1" );

        aParameterLists.set( "bspline_orders", "1,1" );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "initial_refinement", "0,0" );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1;-1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        // aParameterLists.set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        aParameterLists.set( "use_refine_low_level_elements", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );    //,IQIPerimeter_InterfaceVoid,IQIHeatMethodPenalty");
        // aParameterLists.set("output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );
        aParameterLists.set( "PDV_types", "LS1" );

        // initialize fins as swiss cheese geometry
        //( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        // aParameterLists.set( "field_function_name", "Const_Geometry" );
        // aParameterLists.set( "name", "Level_Set_Field" );
        // aParameterLists.set( "number_of_refinements", "0,0" );
        // aParameterLists.set( "refinement_mesh_index", "0,1" );
        // aParameterLists.set( "use_multilinear_interpolation", false );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", -5.0, 2.0, 5.0 );
        aParameterLists.set( "center_y", -1.0, 0.0, 1.0 );
        aParameterLists.set( "normal_x", -2.0, 1.0, 2.0 );
        aParameterLists.set( "normal_y", -1.0, 0.0, 1.0 );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 0.5 );    // FIXME     this has to change

        uint tParamCounter = 0;
        aParameterLists( GEN::PROPERTIES ).add_parameter_list( moris::prm::create_gen_property_parameter_list( gen::Field_Type::CONSTANT ) );
        aParameterLists.set( "name", "LvL_Set_Field" );
        aParameterLists.set( "constant", 0.8 );
        aParameterLists.set( "pdv_type", "LS1" );
        // aParameterLists.set("discretization_mesh_index",   0);
        aParameterLists.set( "discretization_lower_bound", 0.001 );
        aParameterLists.set( "discretization_upper_bound", 1.0 );
        aParameterLists.set( "pdv_mesh_set_names", tTotalDomainAGhost );
        tParamCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tDensityFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerDensityFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "2.0,3.0" );
        aParameterLists.set( "value_function", "tYoungsFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerYoungsFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", std::string( "SPGhost" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tDirichletSets );

        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tVoidInterfaceSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", std::string( "IWGGhost" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists.set( "SOL_save_operator_to_matlab", "Mat.dat" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain );

        aParameterLists.set( "Field_Names", std::string( "UX,UY,VOL" ) );                          //,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT") );
        aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL" ) );                   //,NODAL,NODAL,NODAL,NODAL,NODAL") );
        aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume" ) );    //,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod") );

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
