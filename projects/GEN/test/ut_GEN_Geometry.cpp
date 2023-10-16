/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Geometry.cpp
 *
 */

#include "catch.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_trans.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "fn_check_equal.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"

namespace moris
{

    //------------------------------------------------------------------------------------------------------------------

    uint
    user_defined_phase_function( const ge::Geometry_Bitset& aGeometrySigns )
    {
        uint tPhaseIndex = 3;
        for ( uint tGeometryIndex = 0; tGeometryIndex < 8; tGeometryIndex++ )
        {
            if ( aGeometrySigns.test( tGeometryIndex ) )
            {
                tPhaseIndex = tGeometryIndex;
            }
        }

        return tPhaseIndex;
    }

    //------------------------------------------------------------------------------------------------------------------

    real
    user_defined_geometry_field(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aParameters )
    {
        return aCoordinates( 0 ) * pow( *aParameters( 0 ), 2 ) + aCoordinates( 1 ) * pow( *aParameters( 1 ), 3 );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    user_defined_geometry_sensitivity(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
        aSensitivities = { { 2 * aCoordinates( 0 ) * *aParameters( 0 ), 3 * aCoordinates( 1 ) * pow( *aParameters( 1 ), 2 ) } };
    }

    //------------------------------------------------------------------------------------------------------------------

}

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    // Check for ellipse location in a swiss cheese
    void check_swiss_cheese(
            std::shared_ptr< Level_Set_Geometry > aSwissCheese,
            real                                  aXCenter,
            real                                  aYCenter,
            real                                  aXSemidiameter,
            real                                  aYSemidiameter,
            bool                                  aCheck = true );

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Phase Table", "[gen], [geometry], [phase table]" )
    {
        // Create phase table using number of geometries
        Phase_Table tPhaseTable( 3 );

        // Check number of phases
        CHECK( tPhaseTable.get_num_phases() == 8 );

        // Check individual phases
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 0 ) ) == 0 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 4 ) ) == 1 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 2 ) ) == 2 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 6 ) ) == 3 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 1 ) ) == 4 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 5 ) ) == 5 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 3 ) ) == 6 );
        CHECK( tPhaseTable.get_phase_index( Geometry_Bitset( 7 ) ) == 7 );

        // Create phase custom phase table
        Phase_Table tPhaseTableCustom( 3, Matrix< DDUMat >( { { 3, 2, 1, 0, 0, 1, 2, 3 } } ) );

        // Check number of phases
        CHECK( tPhaseTableCustom.get_num_phases() == 4 );

        // Check individual phases
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 0 ) ) == 3 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 4 ) ) == 2 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 2 ) ) == 1 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 6 ) ) == 0 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 1 ) ) == 0 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 5 ) ) == 1 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 3 ) ) == 2 );
        CHECK( tPhaseTableCustom.get_phase_index( Geometry_Bitset( 7 ) ) == 3 );

        // Create phase custom phase table
        Phase_Table tPhaseTableFunction( &user_defined_phase_function, 4 );

        // Check number of phases
        CHECK( tPhaseTableFunction.get_num_phases() == 4 );

        // Check individual phases
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 0 ) ) == 3 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 4 ) ) == 2 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 2 ) ) == 1 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 6 ) ) == 2 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 1 ) ) == 0 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 5 ) ) == 2 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 3 ) ) == 1 );
        CHECK( tPhaseTableFunction.get_phase_index( Geometry_Bitset( 7 ) ) == 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Circle", "[gen], [geometry], [distributed advs], [circle]" )
    {
        // Set up geometry
        ParameterList tCircle1ParameterList = prm::create_level_set_geometry_parameter_list();
        tCircle1ParameterList.set( "field_type", "circle" );
        tCircle1ParameterList.set( "field_variable_indices", "all" );
        tCircle1ParameterList.set( "adv_indices", "0, 1, 3" );

        ParameterList tCircle2ParameterList = prm::create_level_set_geometry_parameter_list();
        tCircle2ParameterList.set( "field_type", "circle" );
        tCircle2ParameterList.set( "field_variable_indices", "all" );
        tCircle2ParameterList.set( "adv_indices", "0, 2, 4" );

        // ADV vector
        Matrix< DDRMat > tADVs;

        // Distributed ADVs
        sol::Matrix_Vector_Factory tDistributedFactory;

        Matrix< DDSMat >  tADVIds          = { { 0 }, { 1 }, { 2 }, { 3 }, { 4 } };
        sol::Dist_Map*    tADVMap          = tDistributedFactory.create_map( tADVIds );
        sol::Dist_Vector* tDistributedADVs = tDistributedFactory.create_vector( tADVMap, 1, false, true );

        // Define circles
        std::shared_ptr< Level_Set_Geometry > tCircle1;
        std::shared_ptr< Level_Set_Geometry > tCircle2;

        // Loop over ADV types
        for ( bool tDistributed : { false, true } )
        {
            // Set ADVs
            tADVs = { { 0.0, 1.0, 2.0, 1.0, 2.0 } };
            tDistributedADVs->replace_global_values( tADVIds, tADVs );

            // Create circles
            Design_Factory tDesignFactory( { tCircle1ParameterList, tCircle2ParameterList }, tADVs );
            tCircle1 = tDesignFactory.get_geometries()( 0 );
            tCircle2 = tDesignFactory.get_geometries()( 1 );

            // Set distributed ADVs
            if ( tDistributed )
            {
                tCircle1->set_advs( tDistributedADVs );
                tCircle2->set_advs( tDistributedADVs );
            }

            // Set coordinates for checking
            Matrix< DDRMat > tCoordinates0 = { { 0.0, 0.0 } };
            Matrix< DDRMat > tCoordinates1 = { { 1.0, 1.0 } };
            Matrix< DDRMat > tCoordinates2 = { { 2.0, 2.0 } };

            // Check field values
            CHECK( tCircle1->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
            CHECK( tCircle1->get_field_value( 0, tCoordinates1 ) == Approx( 0.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 2.0 ) - 2.0 ) );
            CHECK( tCircle1->get_field_value( 0, tCoordinates2 ) == Approx( sqrt( 5.0 ) - 1.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

            // Check sensitivity values
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { 0.0, 1.0, -1.0 } } ), );
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { 0.0, 1.0, -1.0 } } ), );
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -1.0, 0.0, -1.0 } } ), );
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -sqrt( 2.0 ) / 2.0, sqrt( 2.0 ) / 2.0, -1.0 } } ), );
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -2.0 / sqrt( 5.0 ), -1.0 / sqrt( 5.0 ), -1.0 } } ), );
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -1.0, 0.0, -1.0 } } ), );

            // Change ADVs and coordinates
            tADVs = { { 1.0, 1.0, 2.0, 2.0, 3.0 } };
            tDistributedADVs->replace_global_values( tADVIds, tADVs );
            tCoordinates0( 0 ) = 1.0;
            tCoordinates0( 1 ) = -1.0;
            tCoordinates1( 0 ) = 3.0;
            tCoordinates1( 1 ) = 1.0;
            tCoordinates2( 0 ) = 4.0;
            tCoordinates2( 1 ) = 2.0;

            // Check field values
            CHECK( tCircle1->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
            CHECK( tCircle1->get_field_value( 0, tCoordinates1 ) == Approx( 0.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 5.0 ) - 3.0 ) );
            CHECK( tCircle1->get_field_value( 0, tCoordinates2 ) == Approx( sqrt( 10.0 ) - 2.0 ) );
            CHECK( tCircle2->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

            // Check sensitivity values
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >({ { 0.0, 1.0, -1.0 } }) ,);
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >({ { 0.0, 1.0, -1.0 } }) ,);
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >({ { -1.0, 0.0, -1.0 } }) ,);
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >({ { -2.0 / sqrt( 5.0 ), 1.0 / sqrt( 5.0 ), -1.0 } }) ,);
            CHECK_EQUAL( tCircle1->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >({ { -3.0 / sqrt( 10.0 ), -1.0 / sqrt( 10.0 ), -1.0 } }) ,);
            CHECK_EQUAL( tCircle2->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >({ { -1.0, 0.0, -1.0 } }) ,);
        }

        // Clean up
        delete tDistributedADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Superellipse", "[gen], [geometry], [superellipse]" )
    {
        // Set up geometry
        ParameterList tSuperellipseParameterList = prm::create_level_set_geometry_parameter_list();
        tSuperellipseParameterList.set( "field_type", "superellipse" );
        tSuperellipseParameterList.set( "field_variable_indices", "all" );
        tSuperellipseParameterList.set( "adv_indices", "all" );

        // Create circles
        Matrix< DDRMat > tADVs = { { 3.0, 4.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0 } };
        Design_Factory tDesignFactory( { tSuperellipseParameterList }, tADVs );
        std::shared_ptr< Level_Set_Geometry > tSuperellipse = tDesignFactory.get_geometries()( 0 );

        // Set coordinates for checking
        Matrix< DDRMat > tCoordinates0 = { { 2.0, 2.0 } };
        Matrix< DDRMat > tCoordinates1 = { { 3.0, 3.0 } };
        Matrix< DDRMat > tCoordinates2 = { { 4.0, 4.0 } };

        // Check field values
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates0 ) == Approx( 0.414214 ) );
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates1 ) == Approx( -0.5 ) );
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates0 ),
                Matrix< DDRMat >( { { //
                        7.071067811865476e-01,
                        3.535533905932738e-01,
                        -7.071067811865476e-01,
                        -3.535533905932738e-01,    //
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );

        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates1 ),
                     Matrix< DDRMat >( { { //
                        -0.000000000000000e+00,
                        5.000000000000000e-01,
                        -0.000000000000000e+00,
                        -2.500000000000000e-01,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );

        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates2 ),
                Matrix< DDRMat >( { { //
                        -1.000000000000000e+00,
                        0.000000000000000e+00,
                        -1.000000000000000e+00,
                        -0.000000000000000e+00,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );

        // Change ADVs and coordinates
        tADVs         = { { 2.0, 1.0, 4.0, 3.0, 4.0, 1.0, 0.0, 0.0 } };
        tCoordinates0 = { { -2.0, 1.0 } };
        tCoordinates1 = { { 0.0, 2.5 } };
        tCoordinates2 = { { 2.0, 5.0 } };

        // Check field values
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates1 ) == Approx( pow( 2.0, -0.75 ) - 1.0 ) );
        CHECK( tSuperellipse->get_field_value( 0, tCoordinates2 ) == Approx( 1.0 / 3.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates0 ),
                Matrix< DDRMat >({ { //
                        0.25,
                        0.0,
                        -0.25,
                        0.0,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );

        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates1 ),
                Matrix< DDRMat >({ { //
                        pow( 2.0, 0.25 ) / 8.0,
                        -pow( 2.0, -0.75 ) / 3.0,
                        -pow( 2.0, -0.75 ) / 8.0,
                        -pow( 2.0, -0.75 ) / 6.0,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );

        CHECK_EQUAL( tSuperellipse->get_dfield_dadvs( 0, tCoordinates2 ),
                Matrix< DDRMat >( { { //
                        0.0,
                        -1.0 / 3.0,
                        0.0,
                        -4.0 / 9.0,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX,
                        MORIS_REAL_MAX } } ), );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Sphere", "[gen], [geometry], [sphere]" )
    {
        // Set up geometry
        ParameterList tSphereParameterList = prm::create_level_set_geometry_parameter_list();
        tSphereParameterList.set( "field_type", "sphere" );
        tSphereParameterList.set( "field_variable_indices", "all" );
        tSphereParameterList.set( "adv_indices", "all" );

        // Create sphere
        Matrix< DDRMat > tADVs = { { -1.0, 0.0, 1.0, 2.0 } };
        Design_Factory tDesignFactory( { tSphereParameterList }, tADVs );
        std::shared_ptr< Level_Set_Geometry > tSphere = tDesignFactory.get_geometries()( 0 );

        // Set coordinates for checking
        Matrix< DDRMat > tCoordinates0 = { { 0.0, 0.0, 0.0 } };
        Matrix< DDRMat > tCoordinates1 = { { 1.0, 1.0, 1.0 } };
        Matrix< DDRMat > tCoordinates2 = { { 2.0, 2.0, 2.0 } };

        // Check field values
        CHECK( tSphere->get_field_value( 0, tCoordinates0 ) == Approx( sqrt( 2.0 ) - 2.0 ) );
        CHECK( tSphere->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 5.0 ) - 2.0 ) );
        CHECK( tSphere->get_field_value( 0, tCoordinates2 ) == Approx( sqrt( 14.0 ) - 2.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { -sqrt( 2.0 ) / 2.0, 0.0, sqrt( 2.0 ) / 2.0, -1.0 } } ), );
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -2.0 / sqrt( 5.0 ), -1.0 / sqrt( 5.0 ), 0.0, -1.0 } } ), );
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -3.0 / sqrt( 14.0 ), -sqrt( 2.0 / 7.0 ), -1.0 / sqrt( 14.0 ), -1.0 } } ), );

        // Change ADVs and coordinates
        tADVs         = { { 0.0, 0.0, 1.0, 1.0 } };
        tCoordinates1 = { { 1.0, 1.0, -1.0 } };
        tCoordinates2 = { { 2.0, -2.0, 2.0 } };

        // Check field values
        CHECK( tSphere->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
        CHECK( tSphere->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 6.0 ) - 1.0 ) );
        CHECK( tSphere->get_field_value( 0, tCoordinates2 ) == Approx( 2.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { 0.0, 0.0, 1.0, -1.0 } } ), );
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -1.0 / sqrt( 6.0 ), -1.0 / sqrt( 6.0 ), sqrt( 2.0 / 3.0 ), -1.0 } } ), );
        CHECK_EQUAL( tSphere->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, -1.0 } } ), );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Superellipsoid", "[gen], [geometry], [superellipsoid]" )
    {
        // Set up geometry
        ParameterList tSuperellipsoidParameterList = prm::create_level_set_geometry_parameter_list();
        tSuperellipsoidParameterList.set( "field_type", "superellipsoid" );
        tSuperellipsoidParameterList.set( "field_variable_indices", "all" );
        tSuperellipsoidParameterList.set( "adv_indices", "all" );

        // Create circles
        Matrix< DDRMat > tADVs = { { 3.0, 4.0, 5.0, 1.0, 2.0, 4.0, 3.0 } };
        Design_Factory tDesignFactory( { tSuperellipsoidParameterList }, tADVs );
        std::shared_ptr< Level_Set_Geometry > tSuperellipsoid = tDesignFactory.get_geometries()( 0 );

        // Set coordinates for checking
        Matrix< DDRMat > tCoordinates0 = { { 2.0, 2.0, 5.0 } };
        Matrix< DDRMat > tCoordinates1 = { { 3.0, 3.0, 5.0 } };
        Matrix< DDRMat > tCoordinates2 = { { 4.0, 4.0, 5.0 } };

        // Check field values
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates0 ) == Approx( pow( 2.0, 1.0 / 3.0 ) - 1.0 ) );
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates1 ) == Approx( -0.5 ) );
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates0 ),
                Matrix< DDRMat >( { { //
                        pow( 2.0, -2.0 / 3.0 ),
                        pow( 2.0, -5.0 / 3.0 ),
                        0.0,
                        -pow( 2.0, -2.0 / 3.0 ),
                        -pow( 2.0, -5.0 / 3.0 ),
                        0.0,
                        -0.09703452486 } } ), );

        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates1 ),
                Matrix< DDRMat >( { { //
                        0.0,
                        0.5,
                        0.0,
                        0.0,
                        -0.25,
                        0.0,
                        0.0 } } ), );

        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates2 ),
                Matrix< DDRMat >( { { //
                        -1.0,
                        0.0,
                        0.0,
                        -1.0,
                        0.0,
                        0.0,
                        0.0 } } ), );

        // Change ADVs and coordinates
        tADVs         = { { 2.0, 1.0, 0.0, 5.0, 4.0, 3.0, 4.0 } };
        tCoordinates0 = { { 2.0, -3.0, 0.0 } };
        tCoordinates1 = { { 2.0, -1.0, 1.5 } };
        tCoordinates2 = { { 2.0, 1.0, 4.0 } };

        // Check field values
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates1 ) == Approx( pow( 2.0, -0.75 ) - 1.0 ) );
        CHECK( tSuperellipsoid->get_field_value( 0, tCoordinates2 ) == Approx( 1.0 / 3.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates0 ),
                Matrix< DDRMat >( { { //
                        0.0,
                        0.25,
                        0.0,
                        0.0,
                        -0.25,
                        0.0,
                        0.0 } } ), );

        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates1 ),
                Matrix< DDRMat >( { { //
                        0.0,
                        pow( 2.0, 0.25 ) / 8.0,
                        -pow( 2.0, -0.75 ) / 3.0,
                        0.0,
                        -pow( 2.0, -0.75 ) / 8.0,
                        -pow( 2.0, -0.75 ) / 6.0,
                        -2.575923918612943e-02 } } ), );

        CHECK_EQUAL( tSuperellipsoid->get_dfield_dadvs( 0, tCoordinates2 ),
                Matrix< DDRMat >( { { //
                        0.0,
                        0.0,
                        -1.0 / 3.0,
                        0.0,
                        0.0,
                        -4.0 / 9.0,
                        0.0 } } ), );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "User-defined Geometry", "[gen], [geometry], [user-defined geometry]" )
    {
        // Create user-defined geometry
        Matrix< DDRMat >            tADVs                = { { -1.0, 0.5 } };
        auto tUserDefinedGeometry = std::make_shared< User_Defined_Field >(
                tADVs,
                Matrix< DDUMat >( { { 1, 0 } } ),
                Matrix< DDUMat >( { { 0, 1 } } ),
                Matrix< DDRMat >( { {} } ),
                &user_defined_geometry_field,
                &user_defined_geometry_sensitivity );

        // Set coordinates for checking
        Matrix< DDRMat > tCoordinates1 = { { 1.0, 1.0 } };
        Matrix< DDRMat > tCoordinates2 = { { 2.0, 2.0 } };

        // Check field values
        CHECK( tUserDefinedGeometry->get_field_value( tCoordinates1 ) == Approx( -0.75 ) );
        CHECK( tUserDefinedGeometry->get_field_value( tCoordinates2 ) == Approx( -1.5 ) );

        // Check sensitivity values
        CHECK_EQUAL( tUserDefinedGeometry->get_dfield_dadvs( tCoordinates1 ), Matrix< DDRMat >( { { 1.0, 3.0 } } ), );
        CHECK_EQUAL( tUserDefinedGeometry->get_dfield_dadvs( tCoordinates2 ), Matrix< DDRMat >( { { 2.0, 6.0 } } ), );

        // Change ADVs and coordinates
        tADVs         = { { 2.0, 0.5 } };
        tCoordinates1 = { { 0.0, 1.0 } };
        tCoordinates2 = { { 2.0, -1.0 } };

        // Check field values
        CHECK( tUserDefinedGeometry->get_field_value( tCoordinates1 ) == Approx( 8.0 ) );
        CHECK( tUserDefinedGeometry->get_field_value( tCoordinates2 ) == Approx( -7.5 ) );

        // Check sensitivity values
        CHECK_EQUAL( tUserDefinedGeometry->get_dfield_dadvs( tCoordinates1 ), Matrix< DDRMat >( { { 0.0, 12.0 } } ), );
        CHECK_EQUAL( tUserDefinedGeometry->get_dfield_dadvs( tCoordinates2 ), Matrix< DDRMat >( { { 2.0, -12.0 } } ), );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "B-spline Geometry", "[gen], [geometry], [distributed advs], [B-spline geometry]" )
    {
        // Set up 2 B-spline geometries
        Matrix< DDRMat >                    tADVs( 0, 0 );
        Matrix< DDRMat >                    tRadii = { { 0.5, 0.25 } };
        Cell< std::shared_ptr< Level_Set_Geometry > > tBSplineGeometries( 2 );

        // Loop over possible cases
        for ( uint tBSplineOrder = 1; tBSplineOrder < 3; tBSplineOrder++ )
        {
            for ( uint tLagrangeOrder = 1; tLagrangeOrder < 3; tLagrangeOrder++ )
            {
                for ( uint tRefinement = 0; tRefinement < 3; tRefinement++ )
                {
                    // Create circles
                    for ( uint tGeometryIndex = 0; tGeometryIndex < 2; tGeometryIndex++ )
                    {
                        ParameterList tCircleParameterList = prm::create_level_set_geometry_parameter_list();
                        tCircleParameterList.set( "field_type", "circle" );
                        tCircleParameterList.set( "constant_parameters", "0.0, 0.0, " + std::to_string( tRadii( tGeometryIndex ) ) );
                        tCircleParameterList.set( "discretization_mesh_index", 0 );
                        tCircleParameterList.set( "discretization_lower_bound", -1.0 );
                        tCircleParameterList.set( "discretization_upper_bound", 1.0 );

                        // Set up geometry
                        Design_Factory tDesignFactory( {tCircleParameterList }, tADVs );
                        tBSplineGeometries( tGeometryIndex ) = tDesignFactory.get_geometries()( 0 );
                    }

                    // Create mesh
                    uint tNumElementsPerDimension = 10;
                    mtk::Interpolation_Mesh* tMesh = create_simple_mesh(
                            tNumElementsPerDimension,
                            tNumElementsPerDimension,
                            tLagrangeOrder,
                            tBSplineOrder,
                            tRefinement );

                    // Create geometry engine
                    Geometry_Engine_Parameters tGeometryEngineParameters;
                    tGeometryEngineParameters.mGeometries = tBSplineGeometries;
                    Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                    // Get ADVs and upper/lower bounds
                    tADVs                         = tGeometryEngine.get_advs();
                    Matrix< DDRMat > tLowerBounds = tGeometryEngine.get_lower_bounds();
                    Matrix< DDRMat > tUpperBounds = tGeometryEngine.get_upper_bounds();

                    // Set epsilon for checking
                    real tEpsilon = std::numeric_limits< real >::epsilon() * 10;

                    // Check that ADVs were created and L2 was performed
                    if ( par_rank() == 0 )
                    {
                        uint tNumADVs = pow( tNumElementsPerDimension * pow( 2, tRefinement ) + tBSplineOrder, 2 ) * 2;
                        REQUIRE( tADVs.length() == tNumADVs );
                        REQUIRE( tLowerBounds.length() == tNumADVs );
                        REQUIRE( tUpperBounds.length() == tNumADVs );
                        for ( uint tBSplineIndex = 0; tBSplineIndex < tNumADVs; tBSplineIndex++ )
                        {
                            CHECK( tLowerBounds( tBSplineIndex ) == Approx( -1.0 ) );
                            CHECK( tUpperBounds( tBSplineIndex ) == Approx( 1.0 ) );
                        }
                    }
                    else
                    {
                        REQUIRE( tADVs.length() == 0 );
                        REQUIRE( tLowerBounds.length() == 0 );
                        REQUIRE( tUpperBounds.length() == 0 );
                    }

                    // Epsilon for field value checks must be larger for a quadratic Lagrange mesh
                    if ( tLagrangeOrder > 1 )
                    {
                        tEpsilon = 0.04;
                    }

                    // Quadratic Lagrange, quadratic B-spline is nearly impossible to check
                    if ( tLagrangeOrder < 2 and tBSplineOrder < 2 )
                    {
                        // Looop over both geometries
                        for ( uint tGeometryIndex = 0; tGeometryIndex < 2; tGeometryIndex++ )
                        {
                            // Get geometry back
                            std::shared_ptr< Level_Set_Geometry > tBSplineGeometry = tGeometryEngine.get_geometry( tGeometryIndex );

                            // Get ID Offset
                            uint tOffset = tGeometryIndex * tMesh->get_max_entity_id( mtk::EntityRank::BSPLINE, 0 );

                            // Check field values and sensitivities at all nodes
                            Matrix< DDRMat > tTargetSensitivities;
                            for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
                            {
                                // Get node coordinates
                                Matrix< DDRMat > tNodeCoordinates = tMesh->get_node_coordinate( tNodeIndex );

                                // Set approximate field target
                                Approx tApproxTarget =
                                        Approx( sqrt( pow( tNodeCoordinates( 0 ), 2 ) + pow( tNodeCoordinates( 1 ), 2 ) ) - tRadii( tGeometryIndex ) )
                                                .scale( 2.0 )
                                                .epsilon( tEpsilon );

                                // Check field value
                                CHECK( tBSplineGeometry->get_field_value( tNodeIndex, { {} } ) == tApproxTarget );

                                // Check sensitivities
                                if ( (uint)par_rank() == tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE, 0 ) )
                                {
                                    Matrix< DDRMat > tMatrix = trans( tMesh->get_t_matrix_of_node_loc_ind( tNodeIndex, 0 ) );
                                    Matrix< DDSMat > tIDs    = trans( tMesh->get_coefficient_IDs_of_node( tNodeIndex, 0 ) )
                                                          + tOffset;
                                    CHECK_EQUAL( tBSplineGeometry->get_dfield_dadvs( tNodeIndex, { {} } ), tMatrix ,);
                                    CHECK_EQUAL( tBSplineGeometry->get_determining_adv_ids( tNodeIndex, { {} } ), tIDs ,);
                                }
                            }

                            // Set new ADVs
                            tADVs = tADVs + ( tRadii( tGeometryIndex ) / 2.0 );
                            tGeometryEngine.set_advs( tADVs );

                            // Check field values at all nodes again
                            for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
                            {
                                // Get node coordinates
                                Matrix< DDRMat > tNodeCoordinates = tMesh->get_node_coordinate( tNodeIndex );

                                // Set approximate target
                                Approx tApproxTarget =
                                        Approx( sqrt( pow( tNodeCoordinates( 0 ), 2 ) + pow( tNodeCoordinates( 1 ), 2 ) ) - ( tRadii( tGeometryIndex ) / 2.0 ) )
                                                .scale( 2.0 )
                                                .epsilon( tEpsilon );

                                // Check field value
                                CHECK( tBSplineGeometry->get_field_value( tNodeIndex, { {} } ) == tApproxTarget );
                            }

                            // Reset ADVs for next geometry
                            tADVs = tADVs - ( tRadii( tGeometryIndex ) / 2.0 );
                            tGeometryEngine.set_advs( tADVs );
                        }

                        // Delete mesh pointer
                        delete tMesh;
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Stored Geometry", "[gen], [geometry], [stored geometry]" )
    {
        // Create mesh
        mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 6, 6 );

        // Level set circle parameter list
        ParameterList tCircleParameterList = prm::create_level_set_geometry_parameter_list();
        tCircleParameterList.set( "field_type", "circle" );
        tCircleParameterList.set( "field_variable_indices", "0, 1, 2" );
        tCircleParameterList.set( "adv_indices", "0, 1, 2" );
        tCircleParameterList.set( "discretization_mesh_index", -1 );

        // Set up geometry
        Matrix< DDRMat > tADVs   = { { 0.0, 0.0, 0.5 } };
        Design_Factory tDesignFactory( { tCircleParameterList }, tADVs );
        std::shared_ptr< Level_Set_Geometry > tCircle = tDesignFactory.get_geometries()( 0 );

        // Create geometry engine
        Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = { tCircle };
        Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

        // Get geometry back
        std::shared_ptr< Level_Set_Geometry > tStoredCircle = tGeometryEngine.get_geometry( 0 );

        // Check field values at all nodes
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            // Get node coordinates
            Matrix< DDRMat > tNodeCoordinates = tMesh->get_node_coordinate( tNodeIndex );

            // Check field value
            CHECK( tStoredCircle->get_field_value( tNodeIndex, { {} } ) ==    //
                    Approx( tCircle->get_field_value( tNodeIndex, tNodeCoordinates ) ) );

            // Check sensitivities
            CHECK_EQUAL(
                    tStoredCircle->get_dfield_dadvs( tNodeIndex, { {} } ),
                    tCircle->get_dfield_dadvs( tNodeIndex, tNodeCoordinates ), );
            CHECK_EQUAL(
                    tStoredCircle->get_determining_adv_ids( tNodeIndex, { {} } ),
                    tCircle->get_determining_adv_ids( tNodeIndex, tNodeCoordinates ), );
        }

        // Set new ADVs
        tADVs = { { 1.0, 1.0, 1.0 } };
        tGeometryEngine.set_advs( tADVs );
        tGeometryEngine.reset_mesh_information( tMesh );

        // Check field values at all nodes again
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            // Get node coordinates
            Matrix< DDRMat > tNodeCoordinates = tMesh->get_node_coordinate( tNodeIndex );

            // Check field value
            CHECK( tStoredCircle->get_field_value( tNodeIndex, { {} } ) ==    //
                    Approx( tCircle->get_field_value( tNodeIndex, tNodeCoordinates ) ) );

            // Check sensitivities
            Matrix< DDRMat > tDummyCoordinates = {{}};
            CHECK_EQUAL(
                    tStoredCircle->get_dfield_dadvs( tNodeIndex, tDummyCoordinates ),
                    tCircle->get_dfield_dadvs( tNodeIndex, tNodeCoordinates ), );
            CHECK_EQUAL(
                    tStoredCircle->get_determining_adv_ids( tNodeIndex, tDummyCoordinates ),
                    tCircle->get_determining_adv_ids( tNodeIndex, tNodeCoordinates ), );
        }

        // Delete mesh pointer
        delete tMesh;
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "MultiGeometry", "[gen], [geometry], [combined field]" )
    {
        // ADV indices
        std::string      tADVIndices1    = "0, 1, 3";
        std::string      tADVIndices2    = "0, 2, 4";
        Matrix< DDSMat > tADVIndicesMat1 = string_to_mat< DDSMat >( tADVIndices1 );
        Matrix< DDSMat > tADVIndicesMat2 = string_to_mat< DDSMat >( tADVIndices2 );

        // Set up 2 circles
        Cell< ParameterList > tCircleParameterLists( 2 );
        tCircleParameterLists( 0 ) = prm::create_level_set_geometry_parameter_list();
        tCircleParameterLists( 0 ).set( "field_type", "circle" );
        tCircleParameterLists( 0 ).set( "field_variable_indices", "all" );
        tCircleParameterLists( 0 ).set( "adv_indices", tADVIndices1 );
        tCircleParameterLists( 0 ).set( "name", "circles" );

        tCircleParameterLists( 1 ) = prm::create_level_set_geometry_parameter_list();
        tCircleParameterLists( 1 ).set( "field_type", "circle" );
        tCircleParameterLists( 1 ).set( "field_variable_indices", "all" );
        tCircleParameterLists( 1 ).set( "adv_indices", tADVIndices2 );
        tCircleParameterLists( 1 ).set( "name", "circles" );

        // Create multigeometry
        Matrix< DDRMat > tADVs = { { 0.0, 1.0, 2.0, 1.0, 2.0 } };
        Design_Factory tDesignFactory( tCircleParameterLists, tADVs );
        Cell< std::shared_ptr< Level_Set_Geometry > > tGeometries = tDesignFactory.get_geometries();

        // Should be only one total geometry
        REQUIRE( tGeometries.size() == 1 );
        std::shared_ptr< Level_Set_Geometry > tCombinedField = tGeometries( 0 );

        // Set coordinates for checking
        Matrix< DDRMat > tCoordinates0 = { { 0.0, 0.0 } };
        Matrix< DDRMat > tCoordinates1 = { { 1.0, 1.0 } };
        Matrix< DDRMat > tCoordinates2 = { { 2.0, 2.0 } };

        // Check field values
        CHECK( tCombinedField->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
        CHECK( tCombinedField->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 2.0 ) - 2.0 ) );
        CHECK( tCombinedField->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

        // Check sensitivity values TODO determining IDs
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { 0.0, 1.0, -1.0 } } ) ,);
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -sqrt( 2.0 ) / 2.0, sqrt( 2.0 ) / 2.0, -1.0 } } ) ,);
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -1.0, 0.0, -1.0 } } ) ,);

        // Change ADVs and coordinates
        tADVs( 0 )         = 1.0;
        tADVs( 3 )         = 2.0;
        tADVs( 4 )         = 3.0;
        tCoordinates0( 0 ) = 1.0;
        tCoordinates0( 1 ) = -1.0;
        tCoordinates1( 0 ) = 3.0;
        tCoordinates1( 1 ) = 1.0;
        tCoordinates2( 0 ) = 4.0;
        tCoordinates2( 1 ) = 2.0;

        // Check field values
        CHECK( tCombinedField->get_field_value( 0, tCoordinates0 ) == Approx( 0.0 ) );
        CHECK( tCombinedField->get_field_value( 0, tCoordinates1 ) == Approx( sqrt( 5.0 ) - 3.0 ) );
        CHECK( tCombinedField->get_field_value( 0, tCoordinates2 ) == Approx( 0.0 ) );

        // Check sensitivity values
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates0 ), Matrix< DDRMat >( { { 0.0, 1.0, -1.0 } } ) ,);
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates1 ), Matrix< DDRMat >( { { -2.0 / sqrt( 5.0 ), 1.0 / sqrt( 5.0 ), -1.0 } } ) ,);
        CHECK_EQUAL( tCombinedField->get_dfield_dadvs( 0, tCoordinates2 ), Matrix< DDRMat >( { { -1.0, 0.0, -1.0 } } ) ,);
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Swiss Cheese Slice", "[gen], [geometry], [swiss cheese slice]" )
    {
        SECTION( "Swiss Cheese 1" )
        {
            // Create swiss cheese
            ParameterList tSwissCheeseParameterList = prm::create_swiss_cheese_slice_parameter_list();
            tSwissCheeseParameterList.set( "left_bound", -2.0 );
            tSwissCheeseParameterList.set( "right_bound", 2.0 );
            tSwissCheeseParameterList.set( "bottom_bound", -1.0 );
            tSwissCheeseParameterList.set( "top_bound", 1.0 );
            tSwissCheeseParameterList.set( "hole_x_semidiameter", 0.2 );
            tSwissCheeseParameterList.set( "hole_y_semidiameter", 0.1 );
            tSwissCheeseParameterList.set( "number_of_x_holes", 3 );
            tSwissCheeseParameterList.set( "number_of_y_holes", 5 );
            tSwissCheeseParameterList.set( "superellipse_regularization", 0.0 );
            tSwissCheeseParameterList.set( "superellipse_shift", 0.0 );

            // Create swiss cheese
            Matrix< DDRMat > tADVs = { {} };
            Design_Factory   tDesignFactory( { tSwissCheeseParameterList }, tADVs );
            auto             tSwissCheese = tDesignFactory.get_geometries()( 0 );

            // Check holes
            check_swiss_cheese( tSwissCheese, -2.0, -1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, -2.0, 0.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, -2.0, 1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, -1.0, -1.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, -1.0, 0.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, -1.0, 1.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 0.0, -1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 0.0, 0.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 0.0, 1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 1.0, -1.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 1.0, 0.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 1.0, 1.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 2.0, -1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 2.0, 0.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 2.0, 1.0, 0.2, 0.1 );
            check_swiss_cheese( tSwissCheese, 3.0, -1.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 3.0, 0.0, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.2, 0.1, false );
            check_swiss_cheese( tSwissCheese, 3.0, 1.0, 0.2, 0.1, false );
        }

        SECTION( "Swiss Cheese 2" )
        {
            // Create swiss cheese
            ParameterList tSwissCheeseParameterList = prm::create_swiss_cheese_slice_parameter_list();
            tSwissCheeseParameterList.set( "left_bound", -2.0 );
            tSwissCheeseParameterList.set( "right_bound", 2.0 );
            tSwissCheeseParameterList.set( "bottom_bound", -1.0 );
            tSwissCheeseParameterList.set( "top_bound", 1.0 );
            tSwissCheeseParameterList.set( "number_of_x_holes", 0 );
            tSwissCheeseParameterList.set( "number_of_y_holes", 0 );
            tSwissCheeseParameterList.set( "hole_x_semidiameter", 0.1 );
            tSwissCheeseParameterList.set( "hole_y_semidiameter", 0.2 );
            tSwissCheeseParameterList.set( "target_x_spacing", 2.1 );
            tSwissCheeseParameterList.set( "target_y_spacing", 0.55 );
            tSwissCheeseParameterList.set( "row_offset", 1.0 );

            // Create swiss cheese
            Matrix< DDRMat > tADVs = { {} };
            Design_Factory   tDesignFactory( { tSwissCheeseParameterList }, tADVs );
            auto             tSwissCheese = tDesignFactory.get_geometries()( 0 );

            // Check holes
            check_swiss_cheese( tSwissCheese, -2.0, -1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, -2.0, 0.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, -2.0, 1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, -1.0, -1.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, -1.0, 0.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, -1.0, 1.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 0.0, -1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 0.0, 0.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 0.0, 1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 1.0, -1.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 1.0, 0.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 1.0, 1.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 2.0, -1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 2.0, 0.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 2.0, 1.0, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 3.0, -1.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 3.0, 0.0, 0.1, 0.2, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.1, 0.2 );
            check_swiss_cheese( tSwissCheese, 3.0, 1.0, 0.1, 0.2, false );
        }

        SECTION( "Swiss Cheese 3" )
        {
            // Create swiss cheese
            ParameterList tSwissCheeseParameterList = prm::create_swiss_cheese_slice_parameter_list();
            tSwissCheeseParameterList.set( "left_bound", -2.0 );
            tSwissCheeseParameterList.set( "right_bound", 2.0 );
            tSwissCheeseParameterList.set( "bottom_bound", -1.0 );
            tSwissCheeseParameterList.set( "top_bound", 1.0 );
            tSwissCheeseParameterList.set( "hole_x_semidiameter", 0.45 );
            tSwissCheeseParameterList.set( "hole_y_semidiameter", 0.45 );
            tSwissCheeseParameterList.set( "target_y_spacing", 1.1 );
            tSwissCheeseParameterList.set( "allow_less_than_target_spacing", false );
            tSwissCheeseParameterList.set( "row_offset", 0.0 );

            // Create swiss cheese
            Matrix< DDRMat > tADVs = { {} };
            Design_Factory   tDesignFactory( { tSwissCheeseParameterList }, tADVs );
            auto             tSwissCheese = tDesignFactory.get_geometries()( 0 );

            // Check holes
            check_swiss_cheese( tSwissCheese, -2.0, -1.0, 0.45, 0.45 );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -2.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -2.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -2.0, 1.0, 0.45, 0.45 );
            check_swiss_cheese( tSwissCheese, -1.0, -1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -1.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -1.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, -1.0, 1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 0.0, -1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 0.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 0.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 0.0, 1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 1.0, -1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 1.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 1.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 1.0, 1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 2.0, -1.0, 0.45, 0.45 );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 2.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 2.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 2.0, 1.0, 0.45, 0.45 );
            check_swiss_cheese( tSwissCheese, 3.0, -1.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 3.0, 0.0, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 3.0, -0.5, 0.45, 0.45, false );
            check_swiss_cheese( tSwissCheese, 3.0, 1.0, 0.45, 0.45, false );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    check_swiss_cheese(
            std::shared_ptr< Level_Set_Geometry > aSwissCheese,
            real                        aXCenter,
            real                        aYCenter,
            real                        aXSemidiameter,
            real                        aYSemidiameter,
            bool                        aCheck )
    {
        if ( aCheck )
        {
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter + aXSemidiameter, aYCenter } } ) == Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter, aYCenter + aYSemidiameter } } ) == Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter - aXSemidiameter, aYCenter } } ) == Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter, aYCenter - aYSemidiameter } } ) == Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
        }
        else
        {
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter + aXSemidiameter, aYCenter } } ) != Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter, aYCenter + aYSemidiameter } } ) != Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter - aXSemidiameter, aYCenter } } ) != Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
            CHECK( aSwissCheese->get_field_value( 0, { { aXCenter, aYCenter - aYSemidiameter } } ) != Approx( 0.0 ).margin( 100.0 * MORIS_REAL_EPS ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}
