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
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters )
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
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< real >& aGeometryParameters,
            moris::Matrix< DDRMat >&           aSensitivities )
    {
        // create copy of geometry parameter cell
        moris::Cell< real > tGeometryParameters = aGeometryParameters;

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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

    /* ------------------------------------------------------------------------ */

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 ) / 2.7;

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = ( aCriteria( 1 ) / 1.9 ) - 1.0;

        return tConstraints;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dobjective_dcriteria( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );
        tDObjectiveDCriteria( 0 ) = 1.0 / 2.7;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.numel(), 0.0 );
        tDConstraintDCriteria( 1 ) = 1.0 / 1.9;

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

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", "./Two_Bar_Truss.so" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", "TwoBarTruss.hdf5" );
        tParameterlist( 2 )( 0 ).set( "evaluate_objective_gradients", true );
        tParameterlist( 2 )( 0 ).set( "evaluate_constraint_gradients", true );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "include_bounds", false );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", "1e-6" );
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
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", tInterpolationOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
        //
        //        tParameterlist( 0 )( 0 ).set( "lagrange_input_meshes", "0");
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

        // Main GEN parameter list
        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );

        tParameterlist( 0 )( 0 ).set( "initial_advs",
                std::to_string( 0.0 ) + "," +            //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) );

        tParameterlist( 0 )( 0 ).set( "lower_bounds",
                std::to_string( 0.0 ) + "," +            //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) );

        tParameterlist( 0 )( 0 ).set( "upper_bounds",
                std::to_string( 0.0 ) + "," +            //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 0.0 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.1 ) + "," +    //
                        std::to_string( 1.0 ) + "," +    //
                        std::to_string( 0.5 ) + "," +    //
                        std::to_string( 0.1 ) );

        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy,IQIBulkVolume" );
        tParameterlist( 0 )( 0 ).set( "PDV_types", "" );

        // Geometry parameter lists
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );

        // User defined geometry of bars
        tParameterlist( 1 )( 0 ).set( "field_function_name", "Bars" );
        tParameterlist( 1 )( 0 ).set( "sensitivity_function_name", "BarsGrad" );

        tParameterlist( 1 )( 0 ).set( "field_variable_indices",
                " 0,  1,  2,  3,  4,  5,  6,  7,"
                " 8,  9, 10, 11" );

        tParameterlist( 1 )( 0 ).set( "adv_indices",
                " 0,  1,  2,  3,  4,  5,  6,  7,"
                " 8,  9, 10, 11" );
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

        // properties of bars
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDens );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tEmod );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPois );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties of boundary conditions
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropNeumann" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Neumann" );
        tPropCounter++;

        // properties of bedding (supression for RBMs)
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tBedding );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for outer material
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhost" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBars );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tClampedSurface );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGNeumannFlux" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropNeumann,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLoadedSurface );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG for outer material - ghost
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGGPInnerTemp" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhost,GhostSP" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBarsGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBars );
        tIQICounter++;

        // Max UDisp
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIMaxDofUy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0/2.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBars );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );

        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        tParameterList( 5 )( 0 ).set( "finite_difference_scheme", static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL ) );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", 1e-4 );
    }

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );
        
        tParameterlist( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        //         tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     1 );
        //         tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         1.0 );

        tParameterlist( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        tParameterlist( 7 )( 0 ).set( "ifpack_prec_type", "ILU" );
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
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY,STRAIN_ENERGY,MAXUY,VOLUME" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIBulkStrainEnergy,IQIMaxDofUy,IQIBulkVolume" );
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

