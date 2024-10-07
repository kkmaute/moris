/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Two_Bar_Truss.cpp
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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{

    /* ------------------------------------------------------------------------ */
    // Bar Information
    moris::uint tNumBars  = 2;
    moris::real tElongFac = 1.1;

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    std::string tBars        = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tTotalDomain = tBars;    // + "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tClampedSurface = "SideSet_4_n_p0,SideSet_4_c_p0";
    std::string tLoadedSurface  = "SideSet_2_n_p0,SideSet_2_c_p0";

    std::string tBarsGhost = "ghost_p0";

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string tEmod = "1.0";
    std::string tPois = "0.3";
    std::string tDens = "1.0";

    std::string tBedding = std::to_string( 1.0 * 1e-6 );

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim     = "20, 20";
    std::string tDomainDims         = "1.0, 1.0";
    std::string tDomainOffset       = "0.0, 0.0";
    std::string tDomainSidesets     = "1,2,3,4";
    std::string tInterpolationOrder = "1";

    int tRefineBuffer = 1;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Flag for turning on/off ghost stabilization
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "TwoBarTruss.exo";

    /* ------------------------------------------------------------------------ */

    // Level set function for diamond shaped wedge
    moris::real
    Bars(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        moris::Matrix< DDRMat > tXp = { { aCoordinates( 0 ) }, { aCoordinates( 1 ) }, { 0.0 } };

        // Aux. arrays
        moris::Matrix< moris::DDRMat > tXa( 3, 1 );
        moris::Matrix< moris::DDRMat > tXb( 3, 1 );
        moris::Matrix< moris::DDRMat > tXn( 3, 1 );
        moris::Matrix< moris::DDRMat > tXx( 3, 1 );

        // Variables used to approximate maximum of level set values of bars
        const moris::real tKSbeta = -100.;

        moris::real tKSvalue = 0.0;

        // Loop over all bars
        uint iv = 0;
        for ( uint ib = 0; ib < tNumBars; ++ib )
        {
            // check for correctness of paramter list
            // MORIS_ERROR( (uint) *(aGeometryParameters(iv++)) == ib+1,
            //        "Inconsistency in parameter list: %f  should be %d",*(aGeometryParameters(--iv)),ib);

            tXa( 0 ) = aGeometryParameters( iv++ );
            tXa( 1 ) = aGeometryParameters( iv++ );
            tXa( 2 ) = 0.0;

            moris::real tRada = aGeometryParameters( iv++ );

            tXb( 0 ) = aGeometryParameters( iv++ );
            tXb( 1 ) = aGeometryParameters( iv++ );
            tXb( 2 ) = 0.0;

            moris::real tRadb = aGeometryParameters( iv++ );

            moris::real tLen = norm( tXb - tXa );

            tXn = 1.0 / tLen * ( tXb - tXa );

            moris::Matrix< moris::DDRMat > tSfac = trans( tXn ) * ( tXp - tXa );

            // point along bar
            if ( tSfac( 0 ) >= -tElongFac * tLen && tSfac( 0 ) <= ( tElongFac + 1 ) * tLen )
            {
                moris::real tRad = tRada + ( tRadb - tRada ) * tSfac( 0 ) / tLen;

                moris::real tLscyl = norm( tXp - tXa - tSfac( 0 ) * tXn ) - tRad;

                tKSvalue += std::exp( tKSbeta * tLscyl );
            }
            else
            {
                // point near end with point tXa
                if ( tSfac( 0 ) < -tElongFac * tLen )
                {
                    moris::real tLscap1 = norm( tXp - tXa ) - tRada;

                    tKSvalue += std::exp( tKSbeta * tLscap1 );
                }
                // point near end with point tXb
                else
                {
                    moris::real tLscap2 = norm( tXp - tXb ) - tRadb;

                    tKSvalue += std::exp( tKSbeta * tLscap2 );
                }
            }
        }

        moris::real tVal = 1.0 / tKSbeta * std::log( tKSvalue );

        // clean return value to return non-zero value
        return std::abs( tVal ) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */

    void
    BarsGrad(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters,
            moris::Matrix< DDRMat >&       aSensitivities )
    {
        // create copy of geometry parameter cell
        Vector< real > tGeometryParameters = aGeometryParameters;

        // define finite difference perturbation size
        const real tPerturbation = 1.0e-6;

        // set size of sensitivity vector
        aSensitivities.set_size( 1, tNumBars * 6 );

        // initialize adv counter
        uint iv = 0;

        // loop over all bars
        for ( uint ib = 0; ib < tNumBars; ++ib )
        {
            // loop over all bar advs
            for ( uint ip = 0; ip < 6; ++ip )
            {
                // extract nominal value used for perturbations
                real tNominalValue = tGeometryParameters( iv );

                // set pointer to perturbed value
                tGeometryParameters( iv ) = tNominalValue;

                // positive perturbation
                tNominalValue += tPerturbation;
                real tPosValue = Bars( aCoordinates, tGeometryParameters );

                // positive perturbation
                tNominalValue -= 2.0 * tPerturbation;
                real tNegValue = Bars( aCoordinates, tGeometryParameters );

                // restore nominal value
                tGeometryParameters( iv ) = aGeometryParameters( iv );

                // compute sensitivities
                aSensitivities( iv ) = ( tPosValue - tNegValue ) / ( 2.0 * tPerturbation );

                // increase counter
                iv++;
            }
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Neumann(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // Coordinates of pont
        moris::Matrix< moris::DDRMat > tXp = aFIManager->get_IP_geometry_interpolator()->valx();

        // Check whether point is with loaded surface area
        if ( tXp( 1 ) > 0.475 && tXp( 1 ) < 0.52 )
        {
            aPropMatrix = {
                { 0.0 },
                { 1.0 }
            };
        }
        else
        {
            aPropMatrix = {
                { 0.0 },
                { 0.0 }
            };
        }
    }

    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const( moris::Matrix<
                        moris::DDRMat >&              aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 ) / 2.7;

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = ( aCriteria( 1 ) / 1.9 ) - 1.0;

        return tConstraints;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );
        tDObjectiveDCriteria( 0 ) = 1.0 / 2.7;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );
        tDConstraintDCriteria( 1 ) = 1.0 / 1.9;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", "./Two_Bar_Truss.so" );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        aParameterLists( 2 ).set( "hdf5_path", "TwoBarTruss.hdf5" );
        aParameterLists( 2 ).set( "evaluate_objective_gradients", true );
        aParameterLists( 2 ).set( "evaluate_constraint_gradients", true );
        aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
        aParameterLists( 2 ).set( "include_bounds", false );
        aParameterLists( 2 ).set( "finite_difference_type", "all" );
        aParameterLists( 2 ).set( "finite_difference_epsilons", "1e-6" );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tInterpolationOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists( 0 ).set( "bspline_orders", tInterpolationOrder );
        aParameterLists( 0 ).set( "bspline_pattern", std::string( "0" ) );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "staircase_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
        //
        //        aParameterLists( 0 ).set( "lagrange_input_meshes", "0");
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        // Main GEN parameter list
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        aParameterLists( 0 ).set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );

        // User defined geometry of bars
        aParameterLists( 1 ).set( "field_function_name", "Bars" );
        aParameterLists( 1 ).set( "sensitivity_function_name", "BarsGrad" );
        aParameterLists( 1 )( 0 ).insert( "variable_01", Design_Variable( 0.0, 0.0, 0.0 ) );    // bar1_xa1
        aParameterLists( 1 )( 0 ).insert( "variable_02", Design_Variable( 0.0, 0.0, 0.0 ) );    // bar1_xa2
        aParameterLists( 1 )( 0 ).insert( "variable_03", Design_Variable( 0.1, 0.1, 0.1 ) );    // bar1_ra
        aParameterLists( 1 )( 0 ).insert( "variable_04", Design_Variable( 1.0, 1.0, 1.0 ) );    // bar1_xb1
        aParameterLists( 1 )( 0 ).insert( "variable_05", Design_Variable( 0.5, 0.5, 0.5 ) );    // bar1_xb2
        aParameterLists( 1 )( 0 ).insert( "variable_06", Design_Variable( 0.1, 0.1, 0.1 ) );    // bar1_rb
        aParameterLists( 1 )( 0 ).insert( "variable_07", Design_Variable( 0.0, 0.0, 0.0 ) );    // bar2_xa1
        aParameterLists( 1 )( 0 ).insert( "variable_08", Design_Variable( 1.0, 1.0, 1.0 ) );    // bar2_xa2
        aParameterLists( 1 )( 0 ).insert( "variable_09", Design_Variable( 0.1, 0.1, 0.1 ) );    // bar2_ra
        aParameterLists( 1 )( 0 ).insert( "variable_10", Design_Variable( 1.0, 1.0, 1.0 ) );    // bar2_xb1
        aParameterLists( 1 )( 0 ).insert( "variable_11", Design_Variable( 0.5, 0.5, 0.5 ) );    // bar2_xb2
        aParameterLists( 1 )( 0 ).insert( "variable_12", Design_Variable( 0.1, 0.1, 0.1 ) );    // bar2_rb
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // properties of bars
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDens );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", tEmod );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", tPois );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichlet" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropNeumann" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Neumann" );

        // properties of bedding (supression for RBMs)
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", tBedding );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitsche" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhost" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // create parameter list for IWG 1
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_1" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tBars );

        // create parameter list for IWG 2
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletU" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tClampedSurface );

        // create parameter list for IWG 3
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropNeumann,Traction" );
        aParameterLists( 3 ).set( "mesh_set_names", tLoadedSurface );

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPGhost,GhostSP" );
            aParameterLists( 3 ).set( "mesh_set_names", tBarsGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tBars );

        // Max UDisp
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxDofUy" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "function_parameters", "1.0/2.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tBars );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list

        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
        aParameterLists( 5 ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists( 5 ).set( "finite_difference_perturbation_size", 1e-4 );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        //         aParameterLists( 4 ).set("TSA_Num_Time_Steps",     1 );
        //         aParameterLists( 4 ).set("TSA_Time_Frame",         1.0 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists( 7 ).set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tTotalDomain );
        aParameterLists( 0 ).set( "Field_Names", "UX,UY,STRAIN_ENERGY,MAXUY,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkStrainEnergy,IQIMaxDofUy,IQIBulkVolume" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
