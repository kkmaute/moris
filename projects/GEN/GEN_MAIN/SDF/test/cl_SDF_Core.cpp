/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.cpp
 *
 */

#include <string>
#include <catch.hpp>

// core
#include "typedefs.hpp"
#include "paths.hpp"

// comm
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

// linalg
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "fn_all_true.hpp"

#include "cl_MTK_Mesh_Factory.hpp"

// SDF
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Core.hpp"
#include "cl_SDF_Data.hpp"
#include "cl_SDF_Object.hpp"

using namespace moris;

TEST_CASE(
        "ge::sdf::Generator",
        "[geomeng],[sdf],[Triangle]")
{
    if( par_size() == 1 )
    {
        // get root from environment
        std::string tMorisRoot = moris::get_base_moris_dir();

        // determine path for object file
        std::string tObjectPath = tMorisRoot + "/projects/GEN/GEN_MAIN/SDF/test/data/tetrahedron.obj" ;

        // determine path for mesh file
        std::string tMeshPath =  tMorisRoot + "/projects/GEN/GEN_MAIN/SDF/test/data/TensorMesh.g" ;

        // create triangle object
        sdf::Object tObject( tObjectPath );

        // load MTK mesh from file
        mtk::Mesh * tInput = mtk::create_interpolation_mesh( MeshType::STK, tMeshPath , nullptr);

        // create SDF wrapper for mesh
        sdf::Mesh tMesh( tInput );

        // create data container
        sdf::Data tData( tObject );

        // create core
        sdf::Core tCore( tMesh, tData );

//-------------------------------------------------------------------------------
        SECTION("SDF Core: Raycast Test")
        {
            Matrix< IndexMat > tElementsAtSurfaceExpect;

            if( par_size() == 1 )
            {
                Matrix< IndexMat > tElementsAtSurfaceExpect;
                Matrix< IndexMat > tElementsInVolumeExpect;

                tElementsAtSurfaceExpect =
                {
                    { 73}, { 74}, { 75}, { 76}, { 77}, { 78}, { 81},
                    { 82}, { 83}, { 84}, { 85}, { 86}, { 90}, { 91},
                    { 92}, { 93}, { 98}, { 99}, {100}, {101}, {107},
                    {108}, {115}, {116}, {137}, {138}, {139},
                    {140}, {141}, {142}, {145}, {146}, {149},
                    {150}, {154}, {155}, {156}, {157}, {162},
                    {163}, {164}, {165}, {171}, {172}, {179},
                    {180}, {201}, {202}, {203}, {204}, {205},
                    {206}, {209}, {210}, {211}, {212}, {213},
                    {214}, {218}, {219}, {220}, {221}, {227},
                    {228}, {235}, {236}
                };

                tElementsInVolumeExpect =
                {
                    {147}, {148}
                };

                Matrix< IndexMat > tElementsAtSurface;
                Matrix< IndexMat > tElementsInVolume;

                tCore.calculate_raycast(
                        tElementsAtSurface,
                        tElementsInVolume );

                // print statement needed to put matrix on heap
                print(tElementsInVolumeExpect,"tElementsInVolumeExpect");

                REQUIRE( all_true( tElementsAtSurface == tElementsAtSurfaceExpect ) );
                REQUIRE( all_true( tElementsInVolume == tElementsInVolumeExpect ) );
            }
        }

        // tidy up
        delete tInput;
    }
}

