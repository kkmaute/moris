/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_HMR_User_Defined_Refinement.cpp
 *
 */

#include <catch.hpp>
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src
#include "moris_typedefs.hpp"              // COR/src
#include "cl_Matrix.hpp"                   // LINALG/src
#include "fn_norm.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::hmr
{

    static real LevelSetFunction1( const Matrix< DDRMat >& aPoint )
    {
        return ( norm( aPoint ) - 1.3 ) * ( -1.0 );
    }

    static real LevelSetFunction2( const Matrix< DDRMat >& aPoint )
    {
        return ( norm( aPoint ) - 1.4 ) * ( -1.0 );
    }

    static int user_defined_refinement(
            Element*                aElement,
            const Matrix< DDRMat >& aElementLocalValues )
    {
        int aDoRefine = -1;

        // abs variable field threshold
        real lsth = 0.0;

        // abs variable field bandwidth (absolute)
        real lsbwabs = 0.225;

        // maximum refinement level
        uint maxlevel = 4;

        // min refinement level
        uint minlevel = 0;

        // max refinement level along interface
        uint maxifcref = 4;

        // max refinement level in volume
        uint maxvolref = 1;

        // current refinement level of element
        uint curlevel = aElement->get_level();

        // refinement strategy
        if ( aElementLocalValues.max() >= lsth - lsbwabs )
        {
            // for volume refinement
            if ( aElementLocalValues.min() >= lsth + lsbwabs )
            {
                if ( curlevel < maxvolref && curlevel < maxlevel )
                {
                    aDoRefine = 1;    // refine
                }
                else if ( curlevel == maxvolref || curlevel == minlevel )
                {
                    aDoRefine = 0;    // keep
                }
                else
                {
                    aDoRefine = -1;    // coarsen
                }
            }
            // for interface refinement
            else
            {
                if ( curlevel < maxifcref && curlevel < maxlevel )
                {
                    aDoRefine = 1;    // refine
                }
                else if ( curlevel == maxifcref || curlevel == minlevel )
                {
                    aDoRefine = 0;    // keep
                }
                else
                {
                    aDoRefine = -1;    // coarsen
                }
            }
        }
        else
        {
            if ( curlevel < minlevel )
            {
                aDoRefine = 1;    // refine
            }
            else if ( curlevel == minlevel )
            {
                aDoRefine = 0;    // keep
            }
        }

        return aDoRefine;
    }

    static real Colors(
            const Matrix< DDRMat >& aCoordinates )
    {
        real deps = 1.0e-8;    // offset to avoid machine precision issue when assigning point to geometric object

        real tVal;

        // lower right square
        if ( aCoordinates( 0 ) > 0.1117 + deps && aCoordinates( 1 ) <= -0.11 - deps )
        {
            tVal = 0.0;
        }
        // upper right circle
        else if ( std::sqrt(
                          std::pow( aCoordinates( 0 ) - 2.0, 2 ) + std::pow( aCoordinates( 1 ) - 2.0, 2 ) )
                          - 1.2 - deps
                  <= 0.0 )
        {
            tVal = 1.0;
        }
        // upper left triangle
        else if ( aCoordinates( 0 ) - aCoordinates( 1 ) + 1.5 + deps <= 0.0 )
        {
            tVal = 2.0;
        }
        // lower left circle
        else if ( std::sqrt(
                          std::pow( aCoordinates( 0 ) + 2.0, 2 ) + std::pow( aCoordinates( 1 ) + 2.0, 2 ) )
                          - 2.6 - deps
                  <= 0.0 )
        {
            tVal = 3.0;
        }
        // remainder
        else
        {
            tVal = 4.0;
        }

        return tVal;
    }

    /* ------------------------------------------------------------------------ */

    static int user_defined_refinement_color(
            Element*                aElement,
            const Matrix< DDRMat >& aElementLocalValues )
    {
        int aDoRefine = -1;

        // convert to uint as color is defined here as positive integer
        uint tMaxVal = aElementLocalValues.max();
        uint tMinVal = aElementLocalValues.min();

        if ( tMaxVal != tMinVal )
        {
            aDoRefine = 1;
        }
        else
        {
            aDoRefine = 0;
        }

        return aDoRefine;
    }

    TEST_CASE( "HMR_User_Defined_Refinement", "[moris],[mesh],[hmr],[HMR_User_Defined_Refinement]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 )
        {
            for ( uint tDimension = 2; tDimension <= 2; ++tDimension )
            {
                for ( uint tOrder = 1; tOrder <= 2; tOrder++ )
                {
                    //------------------------------------------------------------------------------
                    //  HMR Parameters setup
                    //------------------------------------------------------------------------------

                    uint tLagrangeMeshIndex = 0;

                    // Dummy parameter list
                    Parameter_List tParam = prm::create_hmr_parameter_list();

                    // The parameter object controls the behavior of HMR.
                    Parameters tParameters;

                    // set values to parameters
                    Matrix< DDLUMat > tNumberOfElements;

                    // set element size
                    tNumberOfElements.set_size( tDimension, 1, 8 );

                    tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

                    tParameters.set_domain_dimensions( { { 2 }, { 2 } } );
                    tParameters.set_domain_offset( { { 0.0 }, { 0.0 } } );

                    // B-Spline truncation is turned on by default.
                    // It is recommended to leave this setting as is.
                    tParameters.set_bspline_truncation( true );

                    tParameters.set_lagrange_orders( { 1 } );
                    tParameters.set_lagrange_patterns( { 0 } );

                    tParameters.set_bspline_orders( { 1 } );
                    tParameters.set_bspline_patterns( { 0 } );

                    tParameters.set_staircase_buffer( 2 );

                    tParameters.set_lagrange_input_mesh( { { 0 } } );

                    tParameters.set_initial_refinement( { { 1 } } );
                    tParameters.set_initial_refinement_patterns( { { 0 } } );

                    tParameters.set_refinement_functions( { &user_defined_refinement } );

                    //------------------------------------------------------------------------------
                    //  HMR Initialization
                    //------------------------------------------------------------------------------

                    // create the HMR object by passing the settings to the constructor
                    HMR tHMR( tParameters );

                    // std::shared_ptr< Database >
                    auto tDatabase = tHMR.get_database();

                    // initial refinement
                    tHMR.perform_initial_refinement();

                    //-----------------------------------------------------------------
                    // First refinement with field 1
                    std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                    // create field
                    std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tLagrangeMeshIndex );

                    tField->evaluate_scalar_function( LevelSetFunction1 );

                    tHMR.based_on_field_put_elements_on_queue( tField->get_node_values(), 0, 0 );

                    tHMR.perform_refinement_based_on_working_pattern( 0, true );

                    //-----------------------------------------------------------------
                    // Second refinement with field 2
                    tField->evaluate_scalar_function( LevelSetFunction2 );

                    tHMR.based_on_field_put_elements_on_queue( tField->get_node_values(), 0, 0 );

                    tHMR.perform_refinement_based_on_working_pattern( 0, true );

                    // calculate T-Matrices etc
                    tDatabase->finalize();

                    tField->evaluate_scalar_function( LevelSetFunction1 );
                    tHMR.save_to_exodus( 0, "UserDefinedRef.exo" );

                    uint tNumElements = tMesh->get_num_elems();

                    // perform test
                    REQUIRE( tNumElements == 1390 );
                }
            }
        }
    }

    TEST_CASE( "Color refinement", "[moris],[mesh],[hmr],[Color refinement]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 )
        {
            uint tDimension = 2;

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            uint tLagrangeMeshIndex = 0;

            // Dummy parameter list
            Parameter_List tParam = prm::create_hmr_parameter_list();

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            // set values to parameters
            Matrix< DDLUMat > tNumberOfElements;

            // set element size
            tNumberOfElements.set_size( tDimension, 1, 20 );

            tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

            tParameters.set_domain_dimensions( { { 4 }, { 4 } } );
            tParameters.set_domain_offset( { { -2.0 }, { -2.0 } } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { 1 } );
            tParameters.set_lagrange_patterns( { 0 } );

            tParameters.set_bspline_orders( { 1 } );
            tParameters.set_bspline_patterns( { 0 } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_initial_refinement( { { 1 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_refinement_functions( { &user_defined_refinement_color } );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // initial refinement
            tHMR.perform_initial_refinement();

            //-----------------------------------------------------------------
            // First refinement with field 1
            std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< Field > tField = tMesh->create_field( "Color", tLagrangeMeshIndex );

            for ( uint Ik = 0; Ik < 3; Ik++ )
            {
                fprintf( stdout, " ----------------------------------------------- \n" );
                tField->evaluate_scalar_function( Colors );

                tHMR.based_on_field_put_elements_on_queue( tField->get_node_values(), 0, 0 );

                tHMR.perform_refinement_based_on_working_pattern( 0, true );
            }

            tField->evaluate_scalar_function( Colors );

            // calculate T-Matrices etc
            tDatabase->finalize();

            tHMR.save_to_exodus( 0, "UserDefinedRefColor.exo" );

            uint tNumElements = tMesh->get_num_elems();

            // perform test
            REQUIRE( tNumElements == 5017 );
        }
    }
}    // namespace moris::hmr
