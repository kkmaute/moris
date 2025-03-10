/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Property.cpp
 *
 */

#include "catch.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "fn_check_equal.hpp"
#include "fn_trans.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Constant property", "[gen], [property], [constant property]" )
    {
        // Set up property
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_gen_property_parameter_list( gen::Field_Type::CONSTANT ) );
        tFieldParameterLists.set( "constant", 0.0, 0.0, 0.0 );
        ADV_Manager tADVManager;
        Design_Factory              tDesignFactory( tFieldParameterLists, tADVManager );
        std::shared_ptr< Property > tConstantProperty = tDesignFactory.get_properties()( 0 );

        // Random distribution
        std::uniform_real_distribution< real > tUniform( -100.0, 100.0 );
        std::default_random_engine             tEngine;

        for ( uint tTestRun = 0; tTestRun < 4; tTestRun++ )
        {
            // Create scaled field
            tADVManager.mADVs( 0 ) = tUniform( tEngine );

            // Loop over coordinate checks
            for ( uint tCoordinateCheck = 0; tCoordinateCheck < 4; tCoordinateCheck++ )
            {
                // Get random coordinates
                Matrix< DDRMat > tCoordinates( { { tUniform( tEngine ), tUniform( tEngine ) } } );

                // Checks
                CHECK( tConstantProperty->get_field_value( 0, tCoordinates ) == Approx( tADVManager.mADVs( 0 ) ) );
                CHECK_EQUAL( tConstantProperty->get_dfield_dadvs( 0, tCoordinates ), { { 1.0 } }, );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Scaled field property", "[gen], [property], [scaled field]" )
    {
        // Set up and create geometry
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tFieldParameterLists.set( "center_x", 0.0, 0.0, 0.0 );
        tFieldParameterLists.set( "center_y", 0.0, 0.0, 0.0 );
        tFieldParameterLists.set( "radius", 0.5, 0.5, 0.5 );
        tFieldParameterLists.set( "name", "My Circle" );

        // Set up property
        tFieldParameterLists.add_parameter_list( prm::create_gen_property_parameter_list( gen::Field_Type::SCALED_FIELD ) );
        tFieldParameterLists.set( "dependencies", "My Circle" );

        // Random distribution
        std::uniform_real_distribution< real > tUniform( -100.0, 100.0 );
        std::default_random_engine             tEngine;

        for ( uint tTestRun = 0; tTestRun < 4; tTestRun++ )
        {
            // Create scaled field
            real tScale = tUniform( tEngine );
            tFieldParameterLists( 1 ).set( "scaling_factor", tScale, false );
            ADV_Manager tADVManager;
            Design_Factory                        tDesignFactory( tFieldParameterLists, tADVManager );
            std::shared_ptr< Level_Set_Geometry > tCircle     = std::dynamic_pointer_cast< Level_Set_Geometry >( tDesignFactory.get_geometries()( 0 ) );
            auto                                  tProperties = tDesignFactory.get_properties();

            // Check that one property was created, and assign it as scaled field
            REQUIRE( tProperties.size() == 1 );
            std::shared_ptr< Property > tScaledField = tProperties( 0 );

            // Loop over coordinate checks
            for ( uint tCoordinateCheck = 0; tCoordinateCheck < 4; tCoordinateCheck++ )
            {
                // Get random coordinates
                Matrix< DDRMat > tCoordinates( { { tUniform( tEngine ), tUniform( tEngine ) } } );

                // Checks
                CHECK( tScaledField->get_field_value( 0, tCoordinates ) == Approx( tCircle->get_field_value( 0, tCoordinates ) * tScale ) );
                CHECK_EQUAL( tScaledField->get_dfield_dadvs( 0, tCoordinates ),
                        Matrix< DDRMat >( tCircle->get_dfield_dadvs( 0, tCoordinates ) * tScale ),
                        1E8, );
                CHECK_EQUAL( tScaledField->get_determining_adv_ids( 0, tCoordinates ),
                        tCircle->get_determining_adv_ids( 0, tCoordinates ), );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "B-spline Property", "[gen], [property], [distributed advs], [B-spline property]" )
    {
        // Constant B-spline parameter list
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_gen_property_parameter_list( gen::Field_Type::CONSTANT ) );
        tFieldParameterLists.set( "constant", 1.0 );
        tFieldParameterLists.set( "discretization_mesh_index", 0 );
        tFieldParameterLists.set( "discretization_lower_bound", -2.0 );
        tFieldParameterLists.set( "discretization_upper_bound", 2.0 );

        // Loop over possible cases
        for ( uint tCaseNumber = 0; tCaseNumber < 4; tCaseNumber++ )
        {
            // Determine mesh orders
            uint tLagrangeOrder = 1;
            uint tBSplineOrder  = 1;
            switch ( tCaseNumber )
            {
                case 1:
                {
                    tLagrangeOrder = 2;
                    break;
                }
                case 2:
                {
                    tBSplineOrder = 2;
                    break;
                }
                case 3:
                {
                    tLagrangeOrder = 2;
                    tBSplineOrder  = 2;
                    break;
                }
                default:
                {
                    // Do nothing
                }
            }

            // Create mesh
            uint                     tNumElementsPerDimension = 10;
            mtk::Interpolation_Mesh* tMesh                    = create_simple_mesh(
                    tNumElementsPerDimension,
                    tNumElementsPerDimension,
                    tLagrangeOrder,
                    tBSplineOrder );

            // Set up property
            ADV_Manager tADVManager;
            Design_Factory   tDesignFactory( tFieldParameterLists, tADVManager );
            auto             tBSplineProperty = tDesignFactory.get_properties()( 0 );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mProperties = { tBSplineProperty };
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Get ADVs and upper/lower bounds
            Vector< real > tADVs = tGeometryEngine.get_advs();
            Vector< real > tLowerBounds = tGeometryEngine.get_lower_bounds();
            Vector< real > tUpperBounds = tGeometryEngine.get_upper_bounds();

            // Check that ADVs were created and L2 was performed
            if ( par_rank() == 0 )
            {
                uint tNumADVs = pow( tNumElementsPerDimension + tBSplineOrder, 2 );
                REQUIRE( tADVs.size() == tNumADVs );
                REQUIRE( tLowerBounds.size() == tNumADVs );
                REQUIRE( tUpperBounds.size() == tNumADVs );
                for ( uint tBSplineIndex = 0; tBSplineIndex < tNumADVs; tBSplineIndex++ )
                {
                    CHECK( tLowerBounds( tBSplineIndex ) == Approx( -2.0 ) );
                    CHECK( tUpperBounds( tBSplineIndex ) == Approx( 2.0 ) );
                }
            }
            else
            {
                REQUIRE( tADVs.size() == 0 );
                REQUIRE( tLowerBounds.size() == 0 );
                REQUIRE( tUpperBounds.size() == 0 );
            }

            // Get property back
            tBSplineProperty = tGeometryEngine.get_property( 0 );

            // Check field values and sensitivities at all nodes
            for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
            {
                // Check field value
                CHECK( tBSplineProperty->get_field_value( tNodeIndex, { {} } ) == Approx( 1.0 ) );

                // Check sensitivities
                if ( (uint)par_rank() == tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE, 0 ) )
                {
                    Matrix< DDRMat > tMatrix = trans( tMesh->get_t_matrix_of_node_loc_ind( tNodeIndex, 0 ) );
                    Matrix< DDSMat > tIDs    = trans( tMesh->get_coefficient_IDs_of_node( tNodeIndex, 0 ) );
                    Vector< sint > tIDVector( tIDs.length() );
                    for ( uint iIndex = 0; iIndex < tIDVector.size(); iIndex++ )
                    {
                        tIDVector( iIndex ) = tIDs( iIndex );
                    }
                    CHECK_EQUAL( tBSplineProperty->get_dfield_dadvs( tNodeIndex, { {} } ), tMatrix, );
                    CHECK_EQUAL( tBSplineProperty->get_determining_adv_ids( tNodeIndex, { {} } ), tIDVector, );
                }
            }

            // Set new ADVs
            for ( uint iADVIndex = 0; iADVIndex < tADVs.size(); iADVIndex++ )
            {
                tADVs( iADVIndex ) = tADVs( iADVIndex ) * 2;
            }
            tGeometryEngine.set_advs( tADVs );

            // Check field values at all nodes again
            for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
            {
                // Check field value
                CHECK( tBSplineProperty->get_field_value( tNodeIndex, { {} } ) == Approx( 2.0 ) );
            }

            // Delete mesh pointer
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
