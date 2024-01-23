/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Surface_Mesh.cpp
 *
 */

#include "catch.hpp"
#include <cmath>
#include "fn_eye.hpp"
#include "paths.hpp"

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_check_equal.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_GEN_Mesh_Field.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Surface Mesh Intersections", "[gen], [pdv], [intersection], [surface mesh intersection]" )
    {
        if ( par_size() == 1 )
        {
            real tEpsilon = 1e-9;

            // create surface mesh geometry
            ParameterList tParameters = prm::create_surface_mesh_geometry_parameter_list();
            tParameters.set( "file_path", moris::get_base_moris_dir() + "projects/GEN/SDF/test/data/rhombus.obj" );
            Surface_Mesh_Parameters                  tSurfaceMeshParameters( tParameters );
            Surface_Mesh_Geometry                    tSurfaceMesh( tSurfaceMeshParameters );
            std::shared_ptr< Surface_Mesh_Geometry > tSurfaceMeshPointer = std::make_shared< Surface_Mesh_Geometry >( tSurfaceMeshParameters );

            // initialize counter for nodes
            moris_index tNodeIndex = 0;

            // create base nodes with parent coordinates
            Matrix< DDRMat > tFirstParentGlobalCoordinates  = { { 0.2, 0.35 } };
            Matrix< DDRMat > tSecondParentGlobalCoordinates = { { 0.3, 0.15 } };

            Background_Node            tFirstBase( tNodeIndex++, tFirstParentGlobalCoordinates );
            Background_Node            tSecondBase( tNodeIndex++, tSecondParentGlobalCoordinates );
            Background_Node            tThirdBase( tNodeIndex++, { { 0.25, 0.15 } } );
            Background_Node            tFourthBase( tNodeIndex++, { { 0.1, 0.2 } } );
            moris::Cell< Node* > tBaseNodes = { &tFirstBase, &tSecondBase, &tThirdBase, &tFourthBase };


            Matrix< DDRMat > tFirstParentParametricCoordinates  = { { -1.0, 1.0 } };
            Matrix< DDRMat > tSecondParentParametricCoordinates = { { 1.0, 1.0 } };

            Background_Node tBaseFirstParent( tNodeIndex++, tFirstParentGlobalCoordinates );
            Background_Node tBaseSecondParent( tNodeIndex++, tSecondParentGlobalCoordinates );

            Parent_Node tFirstParentNode( tBaseFirstParent, tFirstParentParametricCoordinates );
            Parent_Node tSecondParentNode( tBaseSecondParent, tSecondParentParametricCoordinates );

            // create the intersection node
            Intersection_Node* tIntersectionNode = tSurfaceMeshPointer->create_intersection_node(
                    tNodeIndex++,
                    tBaseNodes,
                    tFirstParentNode,
                    tSecondParentNode,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            // check the local coordinate
            real tLocalCoordinateExpected = 0.0;
            real tLocalCoordinate         = tIntersectionNode->get_local_coordinate();
            CHECK( abs( tLocalCoordinateExpected - tLocalCoordinate ) < tEpsilon );
        }
    }

    TEST_CASE( "Engine Surface Mesh Intersections", "[gen], [pdv], [intersection], [surface mesh geometry]," )
    {
        if ( par_size() == 1 )
        {
            // get root from environment
            std::string tMorisRoot = moris::get_base_moris_dir();

            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Set up geometry
            Matrix< DDRMat > tADVs = { {} };

            // surface mesh
            ParameterList tRhombusParameterList = prm::create_surface_mesh_geometry_parameter_list();
            tRhombusParameterList.set( "file_path", tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj" );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mADVs = tADVs;
            Design_Factory tDesignFactory( { tRhombusParameterList }, tADVs );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );
        }
    }
    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::ge
