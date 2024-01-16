/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_SDF_Triangle.cpp
 *
 */

#include <catch.hpp>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "fn_norm.hpp"

#include "cl_SDF_Facet_Vertex.hpp"
#include "cl_SDF_Triangle.hpp"

using namespace moris;
using namespace sdf;

TEST_CASE(
        "ge::sdf::Triangle",
        "[geomeng],[sdf],[Triangle]" )
{
    // example coordinartes for the triangle
    // create list of vertices
    moris::Cell< std::shared_ptr< Facet_Vertex > > tVertices;
    tVertices.resize( 3, nullptr );
    tVertices( 0 ) = std::make_shared< Facet_Vertex >( 0, Matrix< DDRMat >( { { 1.050229216800883 }, { 1.417028287272334 }, { 1.334891429874816 } } ) );
    tVertices( 1 ) = std::make_shared< Facet_Vertex >( 1, Matrix< DDRMat >( { { 0.649117827347835 }, { 1.192051358390624 }, { 0.076610713185436 } } ) );
    tVertices( 2 ) = std::make_shared< Facet_Vertex >( 2, Matrix< DDRMat >( { { 1.088464523227452 }, { -0.179967737969674 }, { 0.309230566913958 } } ) );

    // create single triangle
    Triangle tTriangle( 0, tVertices );

    // epsilon for error
    real tEpsilon  = 1e-9;
    real tEpsilon1 = 1e-6;

    // point needed for some tests
    Matrix< DDRMat > tPoint = {
        { 0.885670154103523 },
        { 0.730289641563979 },
        { 1.121888622244782 }
    };

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: center test" )
    {
        Matrix< F31RMat > tCenter = {
            { 0.929270522458724 },
            { 0.809703969231095 },
            { 0.573577569991403 }
        };

        real tError = norm( tCenter - tTriangle.get_center() );

        REQUIRE( tError < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: normal test" )
    {
        Matrix< DDRMat > tNormal = {
            { -0.912893263456338 },
            { -0.235837187416504 },
            { 0.333176695714916 }
        };
        real tError = norm( tNormal - tTriangle.get_normal() );

        REQUIRE( tError < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Hesse Distance" )
    {
        real tHesse = -0.848180427118634;
        real tError = std::abs( tHesse - tTriangle.get_hesse() );

        REQUIRE( tError < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Bounding Box" )
    {
        // minimum coordinate, expected solution
        Matrix< F31RMat > tMinCoordinate = {
            { 0.649117827347835 },
            { -0.179967737969674 },
            { 0.076610713185436 }
        };

        Matrix< F31RMat > tDelta( 3, 1 );

        for ( uint i = 0; i < 3; ++i )
        {
            tDelta( i ) = tMinCoordinate( i ) - tTriangle.get_min_coord( i );
        }

        REQUIRE( norm( tDelta ) < tEpsilon );

        // maximum coordinate, expected solution
        Matrix< F31RMat > tMaxCoordinate = {
            { 1.088464523227452 },
            { 1.417028287272334 },
            { 1.334891429874816 }
        };

        for ( uint i = 0; i < 3; ++i )
        {
            tDelta( i ) = tMaxCoordinate( i ) - tTriangle.get_max_coord( i );
        }

        REQUIRE( norm( tDelta ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Check Edges" )
    {
        // projection in y-z plane
        REQUIRE( tTriangle.check_edge( 0, 0, tPoint ) == true );
        REQUIRE( tTriangle.check_edge( 1, 0, tPoint ) == false );
        REQUIRE( tTriangle.check_edge( 2, 0, tPoint ) == true );

        // projection in x-z plane
        REQUIRE( tTriangle.check_edge( 0, 1, tPoint ) == true );
        REQUIRE( tTriangle.check_edge( 1, 1, tPoint ) == true );
        REQUIRE( tTriangle.check_edge( 2, 1, tPoint ) == false );

        // projection in x-y plane
        REQUIRE( tTriangle.check_edge( 0, 2, tPoint ) == true );
        REQUIRE( tTriangle.check_edge( 1, 2, tPoint ) == true );
        REQUIRE( tTriangle.check_edge( 2, 2, tPoint ) == true );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Coordinate Axis Intersection Tests" )
    {
        // error that needs to be passed to intersection function
        bool tError;

        // variable for intersection coordinate
        real tXi;

        // expected solution for plane 0
        real tIntersectionYZ = 1.1499023579142249836729065160288;

        // test intersection with plane 0
        tTriangle.intersect_with_coordinate_axis( tPoint, 0, tXi, tError );

        // compare result
        REQUIRE( std::abs( tXi - tIntersectionYZ ) < tEpsilon );

        // expected solution for axis 1
        real tIntersectionXZ = 1.7530961017725150876322698604616;

        // test intersection with axis 1
        tTriangle.intersect_with_coordinate_axis( tPoint, 1, tXi, tError );

        // compare result
        REQUIRE( std::abs( tXi - tIntersectionXZ ) < tEpsilon );

        // expected solution for axis 2
        real tIntersectionXY = 0.39790101461988117245214654257444;

        // test intersection with axis 2
        tTriangle.intersect_with_coordinate_axis( tPoint, 2, tXi, tError );

        REQUIRE( std::abs( tXi - tIntersectionXY ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Area Test" )
    {
        real tArea = 0.974220833569323;

        real tError = std::abs( tArea - tTriangle.get_area() );

        REQUIRE( tError < tEpsilon1 );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Distance Test, part 1" )
    {
        // point coordinates projected into local coordinate system
        Matrix< F31RMat > tProjection = {
            { -0.488600443461426 },
            { 0.109257414478981 },
            { 0.241215798847012 }
        };

        Matrix< F31RMat > tDelta( 3, 1 );
        tDelta = tProjection - tTriangle.project_point_to_local_cartesian( tPoint );

        REQUIRE( norm( tDelta ) < tEpsilon );

        // coordinate of this point in Barycentric coordinates
        Matrix< F31RMat > tBarycentric = {
            { 0.691336629346961 },
            { -0.099792207295336 },
            { 0.408455577948375 }
        };

        tDelta = tBarycentric - tTriangle.get_barycentric_from_local_cartesian( tProjection );

        REQUIRE( norm( tDelta ) < tEpsilon1 );

        // distances from point to each edge
        Matrix< F31RMat > tDistances = {
            { 0.954058443607797 },
            { 0.262060528044059 },
            { 0.641160884691928 }
        };

        for ( uint k = 0; k < 3; ++k )
        {
            tDelta( k ) = tDistances( k )
                        - tTriangle.distance_point_to_edge_in_local_cartesian( tProjection, k );
        }

        REQUIRE( norm( tDelta ) < tEpsilon1 );
    }

    //-------------------------------------------------------------------------------

    SECTION( "SDF Triangle: Distance Test, part 2" )
    {
        // points to test
        Matrix< DDRMat > tPoints = {
            { 0.059234, 0.975223, 0.052644, 0.705246, 0.670342, 0.455441, 1.663529, 0.338271, 1.382485, 1.403031, 1.646998, 1.779594, 0.512393, 0.691528, 0.678431, 1.164209, 0.823881, 1.839212, 0.957501, 1.783941, 1.024442, 0.208111, 1.288813, 1.623248, -0.015943 },
            { -0.113611, 0.088797, 0.975141, 1.195421, 1.682441, 1.012004, 0.694506, 0.804913, 1.378310, 0.413346, 0.208015, 0.130009, 1.783196, 1.433054, 0.773284, 1.513817, 0.104125, 1.447990, 0.393489, 0.810006, 1.232155, 1.255812, 0.866007, 0.752600, 0.015000 },
            { 0.130866, -0.064282, -0.295901, 1.515404, 1.378770, 0.737314, 1.395056, -0.325510, -0.079550, 0.211413, 1.526104, 0.966354, -0.035633, 0.241918, -0.407221, 0.381194, 0.951702, 1.501818, 1.012948, 0.834930, 0.943918, 0.794787, -0.352999, 0.022700, -0.259247 }
        };

        // expected solution
        Matrix< DDRMat > tDistances = {
            { 0.966945530788054 }, { 0.323467402524833 }, { 0.735932134607854 }, { 0.440445580657558 }, { 0.464652493703906 }, { 0.439399071526912 }, { 0.745953762707671 }, { 0.604685436852981 }, { 0.765441776142965 }, { 0.459681594500460 }, { 0.998665297040420 }, { 0.802059606964542 }, { 0.617044825654042 }, { 0.209197658398687 }, { 0.551743281547435 }, { 0.466112521766804 }, { 0.462261824170279 }, { 0.807042092863447 }, { 0.304215233076093 }, { 0.723787800878556 }, { 0.063122262878454 }, { 0.638982157601383 }, { 0.717369656022604 }, { 0.803599693257272 }, { 1.101429307019322 }
        };

        // number of nodes
        uint tNumberOfNodes = tDistances.length();

        Matrix< DDRMat > tError( tNumberOfNodes, 1 );

        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            // create local coordinate from table
            Matrix< F31RMat > tPoint( 3, 1 );
            for ( uint i = 0; i < 3; ++i )
            {
                tPoint( i ) = tPoints( i, k );
            }

            // calculate error
            tError( k ) = tDistances( k ) - tTriangle.get_distance_to_point( tPoint );
        }

        REQUIRE( norm( tError ) < tEpsilon1 );
    }

    //-------------------------------------------------------------------------------
}
