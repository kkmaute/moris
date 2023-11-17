/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HeatMethod.cpp
 *
 */

#include <string>
#include <iostream>
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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// problem dimension: 2D or 3D
extern uint gDim;

// problem dimension: 2D or 3D
extern uint gTestCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    /*
     * FIXME: Set Nitsche penalty to non-zero value for sensitives for cluster
     *        measures have been fixed
     */
    /* ------------------------------------------------------------------------ */

    // 2D or 3D
    bool tIs3D = gDim == 3;

    // turn on/off use of absolute value of level set function
    bool tUseAbsoluteValue = true;

    // interpolation order
    std::string tInterpolationOrder = std::to_string( gInterpolationOrder );

    // prescribed gradient when using non-pdv level set field
    std::string tLevelSetGradxConstant = gDim == 3 ? "1.0;1.0;1.0" : "1.0;1.0";

    // Bspline coefficients for which sensitivities are checked
    std::string tAdvIndicesForFD = gDim == 3 ? "37,61,62,67,93" : "7,8,11,12,13";

    // Nitsche penalty
    std::string tNitschePenalty = "100.0";

    /* ------------------------------------------------------------------------ */

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-7;
    uint tFEMFdScheme  = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );

    // FD step size in sweep
    std::string tFDsweep = "1.0e-7";

    // Number of constraints
    uint tNumConstraints = 6;

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tInnerPhase = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tOuterPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface      = "dbl_iside_p0_1_p1_0";
    std::string tInterfaceInner = "iside_b0_1_b1_0";
    std::string tInterfaceOuter = "iside_b0_0_b1_1";

    std::string tInnerPhaseGhost = "ghost_p1";
    std::string tOuterPhaseGhost = "ghost_p0";

    std::string tTotalDomain = tInnerPhase + "," + tOuterPhase;
    std::string tInterfaces  = tInterface + "," + tInterfaceInner + "," + tInterfaceOuter;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // general
    moris::real tCenterX = 0.0;
    moris::real tCenterY = 0.0;
    moris::real tCenterZ = 0.0;
    moris::real tRadius  = 0.471;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // conductivity
    std::string tConductivity = "1.0";

    // capacity
    std::string tCapacityTheta = "0.1";
    std::string tCapacityPhi   = "0.0";

    // density
    std::string tDensityTheta = "1.0";
    std::string tDensityPhi   = "0.0";

    // prescribed theta on interface
    std::string tPrescTheta = "1.0";

    // prescribed phi on interface
    std::string tPrescPhi = "0.0";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = tIs3D ? "4, 4, 4" : "4, 4";
    std::string tDomainDims     = tIs3D ? "2.0, 2.0, 2.0" : "2.0, 2.0";
    std::string tDomainOffset   = tIs3D ? "-1.0, -1.0, -1.0" : "-1.0, -1.0";
    std::string tDomainSidesets = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";

    std::string tInitialRefinement   = "0";
    std::string tInterfaceRefinement = "0";

    int tRefineBuffer = 0;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 1;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tName          = "HeatMethod";
    std::string tProblemConfig = "_" + std::to_string( gTestCaseIndex );

    std::string tLibraryName    = tName + ".so";
    std::string tOutputFileName = tName + tProblemConfig + ".exo";
    std::string tGENOutputFile  = "GEN_" + tName + tProblemConfig + ".exo";
    std::string tHDF5FileName   = "SEN_" + tName + tProblemConfig + ".hdf5";

    /* ------------------------------------------------------------------------ */

    // Level set function for diamond shaped wedge for GEN
    moris::real
    Inclusion(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
    {
        // distance from sphere center
        moris::real tDx = aCoordinates( 0 ) - tCenterX;
        moris::real tDy = aCoordinates( 1 ) - tCenterY;

        moris::real tVal;
        if ( tIs3D )
        {
            moris::real tDz = aCoordinates( 2 ) - tCenterZ;

            // Compute Signed-Distance field
            tVal = tRadius - std::sqrt( tDx * tDx + tDy * tDy + tDz * tDz );
        }
        else
        {
            // Compute Signed-Distance field
            tVal = tRadius - std::sqrt( tDx * tDx + tDy * tDy );
        }

        // clean return value to return non-zero value
        return std::abs( tVal ) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        // return absolute value of level set function
        aPropMatrix = factor * value;
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of level set function with respect to PDV
    void
    tDerLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->N();
    }

    /* ------------------------------------------------------------------------ */

    // Spatial derivative of level set function defining property in FEM
    void
    tLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->val()( 0 );

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->gradx( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of spatial derivative of level set function with respect to PDV
    void
    tDerLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( ge::PDV_Type::LS1 )->dnNdxn( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const( moris::Matrix<
                        moris::DDRMat >&                   aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( tNumConstraints, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        tObjectives( 0, 0 ) = aCriteria( 0 ) + aCriteria( 1 ) + aCriteria( 2 ) + aCriteria( 3 );

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective      = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "H1Error Const  = " << aCriteria( 0 ) << " ( " << aCriteria( 0 ) / tObjectives( 0, 0 ) << " )\n";
        std::cout << "H1Error PDV    = " << aCriteria( 1 ) << " ( " << aCriteria( 1 ) / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Volume (outer) = " << aCriteria( 2 ) << " ( " << aCriteria( 2 ) / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Volume (outer) = " << aCriteria( 3 ) << " ( " << aCriteria( 3 ) / tObjectives( 0, 0 ) << " )\n";

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );

        for ( uint i = 0; i < tNumConstraints; ++i )
        {
            tConstraints( i ) = aCriteria( i );
        }

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

        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 1.0;
        tDObjectiveDCriteria( 2 ) = 1.0;
        tDObjectiveDCriteria( 3 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.numel(), 0.0 );

        for ( uint i = 0; i < tNumConstraints; ++i )
        {
            tDConstraintDCriteria( i, i ) = 1.0;
        }

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
        tParameterlist( 0 )( 0 ).set( "restart_file", "" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", tHDF5FileName );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tFDsweep );
        tParameterlist( 2 )( 0 ).set( "finite_difference_adv_indices", tAdvIndicesForFD );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", tInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRefinement );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

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
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( moris::prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIH1ErrorConst,IQIH1Error,IQIBulkVolume_Inner,IQIBulkVolume_Outer,IQIThetaMeasure,IQIPhidMeasure" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );

        // Levelset geometry
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        uint tGeoCounter = 0;

        tParameterlist( 1 )( tGeoCounter ).set( "name", "ADVfield" );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Inclusion" );

        tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 0 );
        tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -2.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", 2.0 );
        tGeoCounter++;

        // Levelset property
        tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list() );
        uint tParamCounter = 0;

        tParameterlist( 2 )( tParamCounter ).set( "name", "LevelsetField" );
        tParameterlist( 2 )( tParamCounter ).set( "dependencies", "ADVfield" );
        tParameterlist( 2 )( tParamCounter ).set( "field_type", "scaled_field" );
        tParameterlist( 2 )( tParamCounter ).set( "constant_parameters", "1.0" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomain );
        tParamCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // common properties for theta and phi problems

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tConductivity );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for Theta

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPrescTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPrescTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for phi problem

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPrescPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPrescPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropEigenStrainPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetConst" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetGradxConst" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tLevelSetGradxConstant );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetGradx" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetGradxFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // time continuity weights
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // initial condition
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model - Theta problem
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionTheta" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model - Phi problem
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPhi" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for ghost stabilization parameter for theta and phi problems
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for DBC on interface for theta problem
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", tNitschePenalty );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        //------------------------------------------------------------------------------
        // theta problem

        // create IWG  - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionThetaBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceInner );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceOuter );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerTheta" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhaseGhost );
            tIWGCounter++;

            // create IWG - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPOuterTheta" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhaseGhost );
            tIWGCounter++;
        }

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // theta problem

        // create IWG - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceInner );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceOuter );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerPhi" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInnerPhaseGhost );
            tIWGCounter++;

            // create IWG - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPOuterPhi" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGPTemp,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tOuterPhaseGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal THETA IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTHETA" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Nodal PHID IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkPHID" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Design level set function
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQILevelSet" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Inclusion
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Inner" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInnerPhase );
        tIQICounter++;

        // Matrix
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Outer" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tOuterPhase );
        tIQICounter++;

        // H1 Error if reference is constant
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1ErrorConst" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 0.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // H1 Error if reference is design dependent
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 0.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // integral of theta
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIThetaMeasure" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::MAX_DOF );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0/1.0" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // integral of phid
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIPhidMeasure" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::MAX_DOF );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0/1.0" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();

        tParameterList( 5 )( 0 ).set( "finite_difference_scheme", tFEMFdScheme );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        //------------------------------------------------------------------------------

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        //------------------------------------------------------------------------------

        tParameterlist( 2 ).resize( 3 );
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 2 );

        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 1
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 2 );

        tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 2 );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 4 );
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "THETA" );
        tParameterlist( 3 )( 0 ).set( "NLA_Secondary_DofTypes", "PHID" );

        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 1
        tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "PHID" );
        tParameterlist( 3 )( 1 ).set( "NLA_Secondary_DofTypes", "THETA" );

        tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 2
        tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 1.
        tParameterlist( 3 )( 2 ).set( "NLA_Sub_Nonlinear_Solver", "0,1" );                  // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "THETA;PHID" );

        tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "THETA;PHID" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 2 );    // using NLBGS for forward problem
        // tParameterlist( 4 )( 0 ).set("TSA_nonlinear_solver_for_adjoint_solve", 3 );               // using monlithic for sensitivity problem

        // ----------------------------------------------------------

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "THETA;PHID" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "THETA,1;PHID,1" );

        // ----------------------------------------------------------

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain + "," + tInterfaces );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "THETA,PHID,PHIdesign" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTHETA,IQIBulkPHID,IQILevelSet" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
