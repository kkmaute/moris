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


#include "fn_MTK_mesh_flood_fill.hpp"

// XTKL: Mesh Includes
// #include "cl_MTK_Enums.hpp"

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace moris::mtk
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
        moris_index        tMax = std::numeric_limits< moris_index >::max();
        Matrix< IndexMat > tElementToElement( 7, 3 );
        tElementToElement.fill( tMax );
        tElementToElement( 0, 0 ) = 1;
        tElementToElement( 0, 1 ) = 3;
        tElementToElement( 1, 0 ) = 0;
        tElementToElement( 1, 1 ) = 2;
        tElementToElement( 1, 2 ) = 4;
        tElementToElement( 2, 0 ) = 1;
        tElementToElement( 2, 1 ) = 5;
        tElementToElement( 3, 0 ) = 0;
        tElementToElement( 3, 1 ) = 4;
        tElementToElement( 4, 0 ) = 1;
        tElementToElement( 4, 1 ) = 3;
        tElementToElement( 4, 2 ) = 5;
        tElementToElement( 5, 0 ) = 2;
        tElementToElement( 5, 1 ) = 4;
        tElementToElement( 5, 1 ) = 6;
        tElementToElement( 6, 0 ) = 5;

        // Element Phase Indices
        Vector< moris_index > tElementPhase = { 1, 0, 1, 1, 0, 0, 0 };

        SECTION( "Include all elements" )
        {
            // Active Element Indices (all of them in this case)
            Vector< moris_index > tActiveElements = { 0, 1, 2, 3, 4, 5, 6 };

            Vector< moris_index > tIncludedElementMarker( 7, 1 );    // mark all elements as included

            // Run flood fill Algorithm
            Matrix< IndexMat > tElementSubphase = flood_fill( tElementToElement,
                    tElementPhase,
                    tActiveElements,
                    tIncludedElementMarker,
                    tNumPhases,
                    tMax,
                    true );

            Matrix< IndexMat > tExpElementSubphase( 1, 7 );
            tExpElementSubphase( 0, 0 ) = 0;
            tExpElementSubphase( 0, 1 ) = 1;
            tExpElementSubphase( 0, 2 ) = 2;
            tExpElementSubphase( 0, 3 ) = 0;
            tExpElementSubphase( 0, 4 ) = 1;
            tExpElementSubphase( 0, 5 ) = 1;
            tExpElementSubphase( 0, 6 ) = 1;

            CHECK( all_true( tElementSubphase == tExpElementSubphase ) );
        }
        SECTION( "Excluding element 6" )
        {
            // Active Element Indices (all of them in this case)
            Vector< moris_index > tActiveElements = { 0, 1, 2, 3, 4, 5 };

            Vector< moris_index > tIncludedElementMarker( 7, 1 );    // mark all elements as included
            tIncludedElementMarker( 6 ) = 0;                         // overwrite element 6 and say don't include it

            // Run flood fill Algorithm
            Matrix< IndexMat > tElementSubphase = flood_fill( tElementToElement,
                    tElementPhase,
                    tActiveElements,
                    tIncludedElementMarker,
                    tNumPhases,
                    tMax );

            Matrix< IndexMat > tExpElementSubphase( 1, 6 );
            tExpElementSubphase( 0, 0 ) = 0;
            tExpElementSubphase( 0, 1 ) = 1;
            tExpElementSubphase( 0, 2 ) = 2;
            tExpElementSubphase( 0, 3 ) = 0;
            tExpElementSubphase( 0, 4 ) = 1;
            tExpElementSubphase( 0, 5 ) = 1;

            CHECK( all_true( tElementSubphase == tExpElementSubphase ) );
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
        moris_index tNumPhases = 2;

        // Element to Element Connectivity
        moris_index        tMax = std::numeric_limits< moris_index >::max();
        Matrix< IndexMat > tElementToElement( 7, 3 );
        tElementToElement.fill( tMax );
        tElementToElement( 0, 0 ) = 3;
        tElementToElement( 0, 1 ) = 6;
        tElementToElement( 1, 0 ) = 5;
        tElementToElement( 2, 0 ) = 6;
        tElementToElement( 2, 1 ) = 5;
        tElementToElement( 3, 0 ) = 4;
        tElementToElement( 3, 1 ) = 0;
        tElementToElement( 4, 0 ) = 3;
        tElementToElement( 4, 1 ) = 5;
        tElementToElement( 4, 2 ) = 6;
        tElementToElement( 5, 0 ) = 4;
        tElementToElement( 5, 1 ) = 2;
        tElementToElement( 6, 0 ) = 0;
        tElementToElement( 6, 1 ) = 2;
        tElementToElement( 6, 2 ) = 4;

        // Element Phase Indices
        Vector< moris_index > tElementPhase = { 1, 0, 1, 1, 1, 1, 0 };

        // Active Element Indices (all of them in this case)
        // Note in this case they are not sequential
        Vector< moris_index > tActiveElements = { 0, 6, 2, 3, 4, 5, 1 };

        Vector< moris_index > tIncludedElementMarker( 7, 1 );    // mark all elements as included
        tIncludedElementMarker( 4 ) = 0;                         // overwrite element 4 and say don't include it

        // Run flood fill Algorithm
        Matrix< IndexMat > tElementSubphase = flood_fill( tElementToElement,
                tElementPhase,
                tActiveElements,
                tIncludedElementMarker,
                tNumPhases,
                tMax );

        Matrix< IndexMat > tExpElementSubphase( 1, 6 );
        tExpElementSubphase( 0, 0 ) = 0;
        tExpElementSubphase( 0, 1 ) = 1;
        tExpElementSubphase( 0, 2 ) = 2;
        tExpElementSubphase( 0, 3 ) = 0;
        tExpElementSubphase( 0, 4 ) = 2;
        tExpElementSubphase( 0, 5 ) = 3;

        CHECK( all_true( tElementSubphase == tExpElementSubphase ) );
    }

}    // namespace moris::mtk
