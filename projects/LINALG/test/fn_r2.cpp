/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_r2.cpp
 *
 */

#include <catch.hpp>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_r2.hpp"

TEST_CASE(
        "moris::r2",
        "[linalgebra],[r2]" )
{

    // samples
    moris::Matrix< moris::DDRMat > tFunctionValues = {
            { 0.47168434849678  } ,
            { 0.502683330413384 } ,
            { 0.817697649151299 } ,
            { 0.991287959635897 } ,
            {  0.352454530619332 } ,
            {  0.305813411153125 } ,
            {  -0.238583868483176 } ,
            {  -0.888646523037254 } ,
            {  -0.517977893980116 } ,
            {  -0.024282920100839 } };

    // funciton values
    moris::Matrix< moris::DDRMat > tSamples = {
            { 0.680456229896494 } ,
            {  0.320795927524289 } ,
            {  0.93438246622728 } ,
            {  0.973673990956525 } ,
            {  0.559125096457297 } ,
            {  0.242773944300635 } ,
            { -0.092109813149451 } ,
            { -0.983885560703612 } ,
            { -0.634909924651596 } ,
            { -0.177678956108563 } };

    moris::real tR2 = moris::r2( tFunctionValues, tSamples );

    REQUIRE( std::abs( tR2 - 0.94657872920 ) < 1e-9 );
}

