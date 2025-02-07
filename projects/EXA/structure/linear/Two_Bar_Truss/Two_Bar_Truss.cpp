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

    Vector< uint > tNumElemsPerDim     = { 20, 20 };
    Vector< real > tDomainDims         = { 1.0, 1.0 };
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
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", "./Two_Bar_Truss.so" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", "TwoBarTruss.hdf5" );
        aParameterLists.set( "evaluate_objective_gradients", true );
        aParameterLists.set( "evaluate_constraint_gradients", true );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", "1e-6" );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tInterpolationOrder );
        aParameterLists.set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists.set( "bspline_orders", tInterpolationOrder );
        aParameterLists.set( "bspline_pattern", std::string( "0" ) );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );
        //
        //        aParameterLists.set( "lagrange_input_meshes", "0");
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        // Main GEN parameter list
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );

        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );

        // User defined geometry of bars
        aParameterLists.set( "field_function_name", "Bars" );
        aParameterLists.set( "sensitivity_function_name", "BarsGrad" );
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
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDens );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", tEmod );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", tPois );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropNeumann" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Neumann" );

        // properties of bedding (supression for RBMs)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", tBedding );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        // create parameter list for ghost stabilization parameter for outer material
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // create parameter list for IWG 1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tBars );

        // create parameter list for IWG 2
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tClampedSurface );

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropNeumann,Traction" );
        aParameterLists.set( "mesh_set_names", tLoadedSurface );

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerTemp" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", "SPGhost,GhostSP" );
            aParameterLists.set( "mesh_set_names", tBarsGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists.set( "mesh_set_names", tBars );

        // Max UDisp
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxDofUy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "function_parameters", "1.0/2.0" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tBars );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list

        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", 1e-4 );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
        aParameterLists.set( "preconditioners", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        //         aParameterLists.set("TSA_Num_Time_Steps",     1 );
        //         aParameterLists.set("TSA_Time_Frame",         1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        aParameterLists.set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain );
        aParameterLists.set( "Field_Names", "UX,UY,STRAIN_ENERGY,MAXUY,VOLUME" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkStrainEnergy,IQIMaxDofUy,IQIBulkVolume" );
        aParameterLists.set( "Save_Frequency", 1 );
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
