/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_fn_flood_fill.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "cl_XTK_Enrichment.hpp"

// XTKL: Mesh Includes
#include "cl_MTK_Enums.hpp"

// XTKL: Geometry
#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris::xtk
{

    TEST_CASE( "Generic Floodfill Consecutive Subdomain", "[GEN_FLOOD_FILL]" )
    {
        // This test case colors a 6 element mesh with 2 possible phase indices (0,1)
        // (element id, phase index)
        // *-------*-------*-------*
        // |       |       |       |
        // |  0,1  |  1,0  |  2,1  |
        // |       |       |       |
        // *-------*-------*-------* *-------*
        // |       |       |       | |       |
        // |  3,1  |  4,0  |  5,0  | |  6,0  |
        // |       |       |       | |       |
        // *-------*-------*-------* *-------*
        //
        // Element 6 shows up in the element to element connectivity but I am not interested in considering it
        // This considers a part of the subdomain but the one outside of the subdomain is at the end. Next test will consider
        //
        // Number of phases
        size_t tNumPhases = 2;

        // Element to Element Connectivity
        size_t                          tMax = std::numeric_limits< size_t >::max();
        moris::Matrix< moris::DDSTMat > tElementToElement( 7, 3 );
        tElementToElement.fill( tMax );
        ( tElementToElement )( 0, 0 ) = 1;
        ( tElementToElement )( 0, 1 ) = 3;
        ( tElementToElement )( 1, 0 ) = 0;
        ( tElementToElement )( 1, 1 ) = 2;
        ( tElementToElement )( 1, 2 ) = 4;
        ( tElementToElement )( 2, 0 ) = 1;
        ( tElementToElement )( 2, 1 ) = 5;
        ( tElementToElement )( 3, 0 ) = 0;
        ( tElementToElement )( 3, 1 ) = 4;
        ( tElementToElement )( 4, 0 ) = 1;
        ( tElementToElement )( 4, 1 ) = 3;
        ( tElementToElement )( 4, 2 ) = 5;
        ( tElementToElement )( 5, 0 ) = 2;
        ( tElementToElement )( 5, 1 ) = 4;
        ( tElementToElement )( 5, 1 ) = 6;
        ( tElementToElement )( 6, 0 ) = 5;

        // Element Phase Indices
        moris::Matrix< moris::DDSTMat > tElementPhase( 1, 7 );

        ( tElementPhase )( 0, 0 ) = 1;
        ( tElementPhase )( 0, 1 ) = 0;
        ( tElementPhase )( 0, 2 ) = 1;
        ( tElementPhase )( 0, 3 ) = 1;
        ( tElementPhase )( 0, 4 ) = 0;
        ( tElementPhase )( 0, 5 ) = 0;
        ( tElementPhase )( 0, 6 ) = 0;

        SECTION( "Include all elements" )
        {
            // Active Element Indices (all of them in this case)
            moris::Matrix< moris::DDSTMat > tActiveElements( 1, 7 );
            ( tActiveElements )( 0, 0 ) = 0;
            ( tActiveElements )( 0, 1 ) = 1;
            ( tActiveElements )( 0, 2 ) = 2;
            ( tActiveElements )( 0, 3 ) = 3;
            ( tActiveElements )( 0, 4 ) = 4;
            ( tActiveElements )( 0, 5 ) = 5;
            ( tActiveElements )( 0, 6 ) = 6;

            moris::Matrix< moris::DDSTMat > tIncludedElementMarker( 1, 7 );
            tIncludedElementMarker.fill( 1 );    // Mark all elements as included

            // Run flood fill Algorithm

            moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                    tElementPhase,
                    tActiveElements,
                    tIncludedElementMarker,
                    tNumPhases,
                    tMax,
                    true );

            moris::Matrix< moris::DDSTMat > tExpElementSubphase( 1, 7 );
            ( tExpElementSubphase )( 0, 0 ) = 0;
            ( tExpElementSubphase )( 0, 1 ) = 1;
            ( tExpElementSubphase )( 0, 2 ) = 2;
            ( tExpElementSubphase )( 0, 3 ) = 0;
            ( tExpElementSubphase )( 0, 4 ) = 1;
            ( tExpElementSubphase )( 0, 5 ) = 1;
            ( tExpElementSubphase )( 0, 6 ) = 1;

            CHECK( equal_to( tElementSubphase, tExpElementSubphase ) );
        }
        SECTION( "Excluding element 6" )
        {
            // Active Element Indices (all of them in this case)
            moris::Matrix< moris::DDSTMat > tActiveElements( 1, 6 );
            ( tActiveElements )( 0, 0 ) = 0;
            ( tActiveElements )( 0, 1 ) = 1;
            ( tActiveElements )( 0, 2 ) = 2;
            ( tActiveElements )( 0, 3 ) = 3;
            ( tActiveElements )( 0, 4 ) = 4;
            ( tActiveElements )( 0, 5 ) = 5;

            moris::Matrix< moris::DDSTMat > tIncludedElementMarker( 1, 7 );
            tIncludedElementMarker.fill( 1 );          // Mark all elements as included
            ( tIncludedElementMarker )( 0, 6 ) = 0;    // overwrite element 6 and say don't include it

            // Run flood fill Algorithm
            moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                    tElementPhase,
                    tActiveElements,
                    tIncludedElementMarker,
                    tNumPhases,
                    tMax );

            moris::Matrix< moris::DDSTMat > tExpElementSubphase( 1, 6 );
            ( tExpElementSubphase )( 0, 0 ) = 0;
            ( tExpElementSubphase )( 0, 1 ) = 1;
            ( tExpElementSubphase )( 0, 2 ) = 2;
            ( tExpElementSubphase )( 0, 3 ) = 0;
            ( tExpElementSubphase )( 0, 4 ) = 1;
            ( tExpElementSubphase )( 0, 5 ) = 1;

            CHECK( equal_to( tElementSubphase, tExpElementSubphase ) );
        }
    }

    TEST_CASE( "Generic Floodfill Nonconsectives Subdomain", "[NC_FLOOD_FILL]" )
    {
        // This test case colors a 6 element mesh with 2 possible phase indices (0,1)
        // (element id, phase index)
        // *-------*-------*-------*
        // |       |       |       |
        // |  0,1  |  6,0  |  2,1  |
        // |       |       |       |
        // *-------*-------*-------*-------*
        // |       |  (i)  |       |       |
        // |  3,1  |  4,1  |  5,1  |  1,0  |
        // |       |       |       |       |
        // *-------*-------*-------*-------*
        //
        // Ignoring element 4 (i)

        // Number of phases
        size_t tNumPhases = 2;

        // Element to Element Connectivity
        size_t                          tMax = std::numeric_limits< size_t >::max();
        moris::Matrix< moris::DDSTMat > tElementToElement( 7, 3 );
        tElementToElement.fill( tMax );
        ( tElementToElement )( 0, 0 ) = 3;
        ( tElementToElement )( 0, 1 ) = 6;
        ( tElementToElement )( 1, 0 ) = 5;
        ( tElementToElement )( 2, 0 ) = 6;
        ( tElementToElement )( 2, 1 ) = 5;
        ( tElementToElement )( 3, 0 ) = 4;
        ( tElementToElement )( 3, 1 ) = 0;
        ( tElementToElement )( 4, 0 ) = 3;
        ( tElementToElement )( 4, 1 ) = 5;
        ( tElementToElement )( 4, 2 ) = 6;
        ( tElementToElement )( 5, 0 ) = 4;
        ( tElementToElement )( 5, 1 ) = 2;
        ( tElementToElement )( 6, 0 ) = 0;
        ( tElementToElement )( 6, 1 ) = 2;
        ( tElementToElement )( 6, 2 ) = 4;

        // Element Phase Indices
        moris::Matrix< moris::DDSTMat > tElementPhase( 1, 7 );

        ( tElementPhase )( 0, 0 ) = 1;
        ( tElementPhase )( 0, 1 ) = 0;
        ( tElementPhase )( 0, 2 ) = 1;
        ( tElementPhase )( 0, 3 ) = 1;
        ( tElementPhase )( 0, 4 ) = 1;
        ( tElementPhase )( 0, 5 ) = 1;
        ( tElementPhase )( 0, 6 ) = 0;

        // Active Element Indices (all of them in this case)
        // Note in this case they are not sequential
        moris::Matrix< moris::DDSTMat > tActiveElements( 1, 6 );
        ( tActiveElements )( 0, 0 ) = 0;
        ( tActiveElements )( 0, 1 ) = 6;
        ( tActiveElements )( 0, 2 ) = 2;
        ( tActiveElements )( 0, 3 ) = 3;
        ( tActiveElements )( 0, 4 ) = 5;
        ( tActiveElements )( 0, 5 ) = 1;

        moris::Matrix< moris::DDSTMat > tIncludedElementMarker( 1, 7 );
        tIncludedElementMarker.fill( 1 );          // Mark all elements as included
        ( tIncludedElementMarker )( 0, 4 ) = 0;    // overwrite element 4 and say don't include it

        // Run flood fill Algorithm
        moris::Matrix< moris::DDSTMat > tElementSubphase = flood_fill( tElementToElement,
                tElementPhase,
                tActiveElements,
                tIncludedElementMarker,
                tNumPhases,
                tMax );

        moris::Matrix< moris::DDSTMat > tExpElementSubphase( 1, 6 );
        tExpElementSubphase( 0, 0 ) = 0;
        tExpElementSubphase( 0, 1 ) = 1;
        tExpElementSubphase( 0, 2 ) = 2;
        tExpElementSubphase( 0, 3 ) = 0;
        tExpElementSubphase( 0, 4 ) = 2;
        tExpElementSubphase( 0, 5 ) = 3;

        CHECK( equal_to( tElementSubphase, tExpElementSubphase ) );
    }

}    // namespace moris::xtk
