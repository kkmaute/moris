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
#include "typedefs.hpp"
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            const Cell< real* >&    aGeometryParameters )
    {
        return norm( aCoordinates ) - 3.0;
    }

    //------------------------------------------------------------------------------
    void
    tYoungsFunc( moris::Matrix< moris::DDRMat >&           aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        moris::Matrix< moris::DDRMat > tHCT  = aParameters( 0 );
        moris::real                    tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix = tHCT * std::pow( tDensity, tBeta );
    }

    //------------------------------------------------------------------------------
    void
    tDerYoungsFunc( moris::Matrix< moris::DDRMat >&        aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        moris::Matrix< moris::DDRMat > tHCT  = aParameters( 0 );
        moris::real                    tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        // if( tDensity <0 ){ tDensity = 0.0001; }
        // if( tDensity >1 ){ tDensity = 1; }

        // FIXME density shift missing

        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->N() * tBeta * tHCT( 0 ) * std::pow( tDensity, tBeta - 1 ) / ( 1 - phi_sh );
        // aPropMatrix = tBeta * tHCT * std::pow( tDensity, tBeta -1 ) / ( 1 - phi_sh );
    }

    //------------------------------------------------------------------------------

    void
    tDensityFunc( moris::Matrix< moris::DDRMat >&          aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix.set_size( 1, 1, tDensity );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc( moris::Matrix< moris::DDRMat >&       aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->N() / ( 1 - phi_sh );
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
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
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
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
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
    compute_dobjective_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.numel(), 0.0 );

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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
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
        tParameterlist( 0 )( 0 ).set( "reinitialize_interface_iter", 10000 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", "shape_opt_test.hdf5" );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", "1e-6" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "1,1" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "10,10" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-5.0,-5.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1,1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1,1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0,0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1;-1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        // tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        tParameterlist( 0 )( 0 ).set( "use_refine_low_level_elements", true );
    }

    //--------------------------------------------------------------------------------------------------------------

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
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {

        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy,IQIBulkVolume" );    //,IQIPerimeter_InterfaceVoid,IQIHeatMethodPenalty");
        // tParameterlist(0)(0).set("output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );
        tParameterlist( 0 )( 0 ).set( "PDV_types", "LS1" );
        tParameterlist( 0 )( 0 ).set( "initial_advs", "2.0; 0.0; 1.0; 0.0" );
        tParameterlist( 0 )( 0 ).set( "lower_bounds", "-5.0; -1.0; -2.0; -1.0" );
        tParameterlist( 0 )( 0 ).set( "upper_bounds", "5.0; 1.0; 2.0; 1.0" );

        // init geometry counter
        uint tGeoCounter = 0;

        // initialize fins as swiss cheese geometry
        //( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        // tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Const_Geometry" );
        // tParameterlist( 1 )( tGeoCounter ).set( "name", "Level_Set_Field" );
        // tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "0,0" );
        // tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        // tParameterlist( 1 )( tGeoCounter ).set( "multilinear_intersections", false );

        tParameterlist( 1 ).push_back( prm::create_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "type", "plane" );
        // tParameterlist(1)( tGeoCounter ).set("constant_parameters", "0.0, 0.0, 2.0");
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 10e-14 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.5 );    // FIXME     this has to change
        tParameterlist( 1 )( tGeoCounter ).set( "field_variable_indices", "0, 1, 2, 3" );
        tParameterlist( 1 )( tGeoCounter ).set( "adv_indices", "0, 1, 2, 3" );

        uint tParamCounter = 0;
        tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list() );
        tParameterlist( 2 )( tParamCounter ).set( "name", "LvL_Set_Field" );
        tParameterlist( 2 )( tParamCounter ).set( "type", "constant" );
        tParameterlist( 2 )( tParamCounter ).set( "constant_parameters", "0.8" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
        // tParameterlist( 2 )( tParamCounter ).set("discretization_mesh_index",   0);
        tParameterlist( 2 )( tParamCounter ).set( "discretization_lower_bound", 0.001 );
        tParameterlist( 2 )( tParamCounter ).set( "discretization_upper_bound", 1.0 );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomainAGhost );
        tParamCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
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
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tDensityFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerDensityFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0,3.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tYoungsFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerYoungsFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
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
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulkSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tDirichletSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVoidInterfaceSets );
        tIWGCounter++;

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhost" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 ).resize( 1 );
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 5 );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 1 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );
        tParameterlist( 3 )( 0 ).set( "NLA_Secondary_DofTypes", "" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 0 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "Mat.dat" );

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "UX", 0 );
        tParameterlist( 0 )( 0 ).set( "UY", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain );

        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,VOL" ) );                          //,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT") );
        tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL" ) );                   //,NODAL,NODAL,NODAL,NODAL,NODAL") );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume" ) );    //,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod") );

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
