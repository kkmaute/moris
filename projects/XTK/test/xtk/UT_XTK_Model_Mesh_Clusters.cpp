/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Model_Mesh_Clusters.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Sphere.hpp"

namespace moris::xtk
{

    TEST_CASE( "Mesh Cluster Output", "[XTK] [XTK_CLUSTER]" )
    {
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
        MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

        if ( tProcSize <= 4 )
        {
            // Geometry Engine Setup ---------------------------------------------------------
            // Using a Levelset Sphere as the Geometry

            real                                              tRadius   = 0.25;
            real                                              tXCenter  = 1.0;
            real                                              tYCenter  = 1.0;
            real                                              tZCenter  = 0.0;
            auto                                              tSphere   = std::make_shared< moris::gen::Sphere >( tXCenter, tYCenter, tZCenter, tRadius );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tSphere ) };

            // Create Mesh --------------------------------------------------------------------
            std::string                     tMeshFileName = "generated:1x1x4|sideset:Z";
            moris::mtk::Interpolation_Mesh* tMeshData     = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, NULL );

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model  tXTKModel( tModelDimension, tMeshData, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose( tDecompositionMethods );

            delete tMeshData;
        }
    }
}    // namespace moris::xtk
