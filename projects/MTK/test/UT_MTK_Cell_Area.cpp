/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Cell_Area.cpp
 *
 */

#include "catch.hpp"
#include "cl_Communication_Tools.hpp"

// base class
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Cell_STK.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Vertex_STK.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Space_Interpolator.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_reindex_mat.hpp"

namespace moris::mtk
{
    TEST_CASE( "MTK Cell Size Integration", "[MTK],[Size_Integration],[Size_Quad4]" )
    {
        if ( par_size() <= 1 )
        {
            // create a 2D MORIS mesh of quad4's using MTK database
            //------------------------------------------------------------------------------
            uint tNumDim = 2;    // specify number of spatial dimensions

            // Node coordinate matrix
            Matrix< DDRMat > tCoords = { { 0.0, 0.0 },
                { 1.0, 0.0 },
                { 1.7, 1.7 },
                { 0.7, 1.5 } };

            Matrix< IndexMat > tNodeIndices = { { 0, 1, 2, 3 } };
            Matrix< IndexMat > tNodeIds     = { { 1, 2, 3, 4 } };

            // specify element connectivity of quad for mesh
            Matrix< IdMat > tElementConnQuad = { { 1, 2, 3, 4 } };

            // specify the local to global element map for quads
            Matrix< IdMat > tElemLocalToGlobalQuad = { { 1 } };

            // specify the local to global map
            // Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &tNumDim;
            aMeshData.ElemConn( 0 )             = &tElementConnQuad;
            aMeshData.CellTopology( 0 )         = CellTopology::QUAD4;
            aMeshData.NodeCoords                = &tCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocalToGlobalQuad;
            // aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

            Mesh *tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            // Dump to file
            std::string tFileOutput = "./mtk_quad4_cell_area_ut.exo";
            tMesh->create_output_mesh( tFileOutput );

            CHECK( tMesh->get_num_elems() == 1 );

            // Setup Node Vertices (note: this data structure will be in the STK_Implementation
            Vector< Vertex * > tElementVertices;
            for ( size_t i = 0; i < tNodeIndices.numel(); i++ )
            {
                tElementVertices.push_back( new Vertex_STK( tNodeIds( i ), tNodeIndices( i ), tMesh ) );
            }

            // Setup cell associated with element index 0
            std::shared_ptr< Cell_Info > tQuad4 = std::make_shared< Cell_Info_Quad4 >();
            Cell_STK                     tCell( tQuad4, 1, 0, tElementVertices, tMesh );

            // Some checks on the cell
            CHECK( tCell.get_id() == 1 );
            CHECK( tCell.get_index() == 0 );
            CHECK( tCell.get_owner() == (moris_id)par_rank() );
            CHECK( tCell.get_number_of_vertices() == 4 );
            Interpolation_Order tInterpOrder = tCell.get_interpolation_order();
            CHECK( tInterpOrder == Interpolation_Order::LINEAR );

            Integration_Order tIntegOrder = tCell.get_integration_order();
            CHECK( tIntegOrder == Integration_Order::QUAD_2x2 );

            CHECK( tCell.get_cell_info()->compute_cell_size_general( &tCell ) == Approx( 1.53 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_general( &tCell ) == tCell.get_cell_info()->compute_cell_size_straight( &tCell ) );
        }
    }

    TEST_CASE( "MTK Cell Size Integration Hex", "[MTK],[Size_Integration],[Size_Hex8]" )
    {
        if ( par_size() <= 1 )
        {
            // create a 2D MORIS mesh of quad4's using MTK database
            //------------------------------------------------------------------------------
            uint tNumDim = 3;    // specify number of spatial dimensions

            // Node coordinate matrix
            Matrix< DDRMat > tCoords            = { { 0.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.5 },
                           { 1.0, 1.0, 0.5 },
                           { 0.0, 1.0, 0.0 },
                           { 0.0, 0.3, 1.0 },
                           { 1.0, 0.3, 1.5 },
                           { 1.0, 1.3, 1.5 },
                           { 0.0, 1.3, 1.0 } };
            Matrix< IdMat >  tNodeLocalToGlobal = { { 1, 2, 3, 4, 5, 6, 7, 8 } };

            // specify element connectivity of quad for mesh
            Matrix< IdMat > tElementConnQuad = { { 1, 2, 3, 4, 5, 6, 7, 8 } };

            // specify the local to global element map for quads
            Matrix< IdMat > tElemLocalToGlobalQuad = { { 1 } };

            Matrix< IndexMat > tNodeIndices = { { 0, 1, 2, 3, 4, 5, 6, 7 } };
            Matrix< IndexMat > tNodeIds     = { { 1, 2, 3, 4, 5, 6, 7, 8 } };

            // specify the local to global map
            // Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &tNumDim;
            aMeshData.ElemConn( 0 )             = &tElementConnQuad;
            aMeshData.CellTopology( 0 )         = CellTopology::HEX8;
            aMeshData.NodeCoords                = &tCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocalToGlobalQuad;
            aMeshData.LocaltoGlobalNodeMap      = &tNodeLocalToGlobal;

            Mesh *tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            // Dump to file
            std::string tFileOutput = "./mtk_hex8_cell_area_ut.exo";
            tMesh->create_output_mesh( tFileOutput );

            CHECK( tMesh->get_num_elems() == 1 );

            // Setup Node Vertices (note: this data structure will be in the STK_Implementation
            Vector< Vertex * > tElementVertices;
            for ( size_t i = 0; i < tNodeIndices.numel(); i++ )
            {
                tElementVertices.push_back( new Vertex_STK( tNodeIds( i ), tNodeIndices( i ), tMesh ) );
            }

            // Setup cell associated with element index 0
            std::shared_ptr< Cell_Info > tHex8 = std::make_shared< Cell_Info_Hex8 >();
            Cell_STK                     tCell( tHex8, 1, 0, tElementVertices, tMesh );

            // Some checks on the cell
            CHECK( tCell.get_id() == 1 );
            CHECK( tCell.get_index() == 0 );
            CHECK( tCell.get_owner() == (moris_id)par_rank() );
            CHECK( tCell.get_number_of_vertices() == 8 );
            CHECK( tCell.get_cell_info()->compute_cell_size_general( &tCell ) == Approx( 1.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_straight( &tCell ) == Approx( 1.0 ) );
        }
    }

    TEST_CASE( "MTK Quad4 Size Derivative", "[MTK],[Size_Deriv],[Size_Deriv_Quad4]" )
    {
        if ( par_size() <= 1 )
        {
            // create a 2D MORIS mesh of quad4's using MTK database
            //------------------------------------------------------------------------------
            uint tNumDim = 2;    // specify number of spatial dimensions

            // Node coordinate matrix
            Matrix< DDRMat > tCoords = { { 0.0, 0.0 },
                { 1.0, 0.0 },
                { 1.0, 1.0 },
                { 0.0, 1.0 } };

            Matrix< IndexMat > tNodeIndices = { { 0, 1, 2, 3 } };
            Matrix< IndexMat > tNodeIds     = { { 1, 2, 3, 4 } };

            // specify element connectivity of quad for mesh
            Matrix< IdMat > tElementConnQuad = { { 1, 2, 3, 4 } };

            // specify the local to global element map for quads
            Matrix< IdMat > tElemLocalToGlobalQuad = { { 1 } };

            // specify the local to global map
            // Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &tNumDim;
            aMeshData.ElemConn( 0 )             = &tElementConnQuad;
            aMeshData.CellTopology( 0 )         = CellTopology::QUAD4;
            aMeshData.NodeCoords                = &tCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocalToGlobalQuad;
            // aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

            Mesh *tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );

            CHECK( tMesh->get_num_elems() == 1 );

            // Setup Node Vertices (note: this data structure will be in the STK_Implementation
            Vector< Vertex * > tElementVertices;
            for ( size_t i = 0; i < tNodeIndices.numel(); i++ )
            {
                tElementVertices.push_back( new Vertex_STK( tNodeIds( i ), tNodeIndices( i ), tMesh ) );
            }

            // Setup cell associated with element index 0
            std::shared_ptr< Cell_Info > tQuad4 = std::make_shared< Cell_Info_Quad4 >();
            Cell_STK                     tCell( tQuad4, 1, 0, tElementVertices, tMesh );

            // Some checks on the cell
            CHECK( tCell.get_id() == 1 );
            CHECK( tCell.get_index() == 0 );
            CHECK( tCell.get_owner() == (moris_id)par_rank() );
            CHECK( tCell.get_number_of_vertices() == 4 );
            Interpolation_Order tInterpOrder = tCell.get_interpolation_order();
            CHECK( tInterpOrder == Interpolation_Order::LINEAR );

            Integration_Order tIntegOrder = tCell.get_integration_order();
            CHECK( tIntegOrder == Integration_Order::QUAD_2x2 );

            // checking derivative of 0 vertex in the 0 direction
            CHECK( tCell.compute_cell_measure_deriv( 0, 0 ) == Approx( -0.5 ) );

            // checking derivative of 0 vertex in the 1 direction
            CHECK( tCell.compute_cell_measure_deriv( 0, 1 ) == Approx( -0.5 ) );

            // checking other nodes
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 1 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 1 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 0 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 1 ) == Approx( 0.5 ) );

            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 0 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 1 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 1 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 1 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 0 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 1 ) == Approx( 0.5 ) );
        }
    }

    TEST_CASE( "MTK Tri3 Size Derivative", "[MTK],[Size_Deriv],[Size_Deriv_Tri3]" )
    {
        if ( par_size() <= 1 )
        {
            // create a 2D MORIS mesh of quad4's using MTK database
            //------------------------------------------------------------------------------
            uint tNumDim = 2;    // specify number of spatial dimensions

            // Node coordinate matrix
            Matrix< DDRMat > tCoords = { { 0.0, 0.0 },
                { 1.0, 0.0 },
                { 0.0, 1.0 } };

            Matrix< IndexMat > tNodeIndices = { { 0, 1, 2 } };
            Matrix< IndexMat > tNodeIds     = { { 1, 2, 3 } };

            // specify element connectivity of quad for mesh
            Matrix< IdMat > tElementConnQuad = { { 1, 2, 3 } };

            // specify the local to global element map for quads
            Matrix< IdMat > tElemLocalToGlobalQuad = { { 1 } };

            // specify the local to global map
            // Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &tNumDim;
            aMeshData.ElemConn( 0 )             = &tElementConnQuad;
            aMeshData.CellTopology( 0 )         = CellTopology::TRI3;
            aMeshData.NodeCoords                = &tCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocalToGlobalQuad;
            // aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

            Mesh *tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );

            CHECK( tMesh->get_num_elems() == 1 );

            // Setup Node Vertices (note: this data structure will be in the STK_Implementation
            Vector< Vertex * > tElementVertices;
            for ( size_t i = 0; i < tNodeIndices.numel(); i++ )
            {
                tElementVertices.push_back( new Vertex_STK( tNodeIds( i ), tNodeIndices( i ), tMesh ) );
            }

            // Setup cell associated with element index 0

            std::shared_ptr< Cell_Info > tTri3 = std::make_shared< Cell_Info_Tri3 >();

            Cell_STK tCell( tTri3, 1, 0, tElementVertices, tMesh );

            // Some checks on the cell
            CHECK( tCell.get_id() == 1 );
            CHECK( tCell.get_index() == 0 );
            CHECK( tCell.get_owner() == (moris_id)par_rank() );
            CHECK( tCell.get_number_of_vertices() == 3 );
            Interpolation_Order tInterpOrder = tCell.get_interpolation_order();
            CHECK( tInterpOrder == Interpolation_Order::LINEAR );

            // checking derivative of 0 vertex in the 0 direction
            CHECK( tCell.compute_cell_measure_deriv( 0, 0 ) == Approx( -0.5 ) );

            // checking derivative of 0 vertex in the 1 direction
            CHECK( tCell.compute_cell_measure_deriv( 0, 1 ) == Approx( -0.5 ) );

            // checking other nodes
            CHECK( tCell.compute_cell_measure_deriv( 1, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.compute_cell_measure_deriv( 1, 1 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 2, 0 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 2, 1 ) == Approx( 0.5 ) );

            // using general calc
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 0 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 1 ) == Approx( -0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 0 ) == Approx( 0.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 1 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 0 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 1 ) == Approx( 0.5 ) );
        }
    }

    TEST_CASE( "MTK Tet4 Size Derivative", "[MTK],[Size_Deriv],[Size_Deriv_Tet4]" )
    {
        if ( par_size() <= 1 )
        {
            // create a 2D MORIS mesh of tet4's using MTK database
            //------------------------------------------------------------------------------
            uint tNumDim = 3;    // specify number of spatial dimensions

            // Node coordinate matrix
            Matrix< DDRMat > tCoords = { { 0.0, 0.0, 0.0 },
                { 3.0, 0.0, 0.0 },
                { 0.0, 3.0, 0.0 },
                { 0.0, 0.0, 6.0 } };

            Matrix< IndexMat > tNodeIndices = { { 0, 1, 2, 3 } };
            Matrix< IndexMat > tNodeIds     = { { 1, 2, 3, 4 } };

            // specify element connectivity of tet for mesh
            Matrix< IdMat > tElementConnQuad = { { 1, 2, 3, 4 } };

            // specify the local to global element map for tets
            Matrix< IdMat > tElemLocalToGlobalQuad = { { 1 } };

            // specify the local to global map
            // Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &tNumDim;
            aMeshData.ElemConn( 0 )             = &tElementConnQuad;
            aMeshData.CellTopology( 0 )         = CellTopology::TET4;
            aMeshData.NodeCoords                = &tCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocalToGlobalQuad;
            // aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

            Mesh *tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );

            CHECK( tMesh->get_num_elems() == 1 );

            // Setup Node Vertices (note: this data structure will be in the STK_Implementation
            Vector< Vertex * > tElementVertices;
            for ( size_t i = 0; i < tNodeIndices.numel(); i++ )
            {
                tElementVertices.push_back( new Vertex_STK( tNodeIds( i ), tNodeIndices( i ), tMesh ) );
            }

            // Setup cell associated with element index 0
            std::shared_ptr< Cell_Info > tTet4 = std::make_shared< Cell_Info_Tet4 >();

            Cell_STK tCell( tTet4, 1, 0, tElementVertices, tMesh );

            // Some checks on the cell
            CHECK( tCell.get_id() == 1 );
            CHECK( tCell.get_index() == 0 );
            CHECK( tCell.get_owner() == (moris_id)par_rank() );
            CHECK( tCell.get_number_of_vertices() == 4 );
            Interpolation_Order tInterpOrder = tCell.get_interpolation_order();
            CHECK( tInterpOrder == Interpolation_Order::LINEAR );

            // check the area of the cell
            CHECK( tCell.get_cell_info()->compute_cell_size_straight( &tCell ) == Approx( 9.0 ) );
            CHECK( tCell.compute_cell_measure() == Approx( 9.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_general( &tCell ) == Approx( 9.0 ) );

            // checking derivative in the 0 direction of nodes 0-3
            CHECK( tCell.compute_cell_measure_deriv( 0, 0 ) == Approx( -3.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 1, 0 ) == Approx( 3.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 2, 0 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 3, 0 ) == Approx( 0.0 ) );

            // checking derivative in the 1 direction of nodes 0-3
            CHECK( tCell.compute_cell_measure_deriv( 0, 1 ) == Approx( -3.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 1, 1 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 2, 1 ) == Approx( 3.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 3, 1 ) == Approx( 0.0 ) );

            // checking derivative in the 2 direction of nodes 0-3
            CHECK( tCell.compute_cell_measure_deriv( 0, 2 ) == Approx( -1.5 ) );
            CHECK( tCell.compute_cell_measure_deriv( 1, 2 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 2, 2 ) == Approx( 0.0 ) );
            CHECK( tCell.compute_cell_measure_deriv( 3, 2 ) == Approx( 1.5 ) );

            // checking derivative in the 0 direction of nodes 0-3
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 0 ) == Approx( -3.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 0 ) == Approx( 3.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 0 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 0 ) == Approx( 0.0 ) );

            // checking derivative in the 1 direction of nodes 0-3
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 1 ) == Approx( -3.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 1 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 1 ) == Approx( 3.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 1 ) == Approx( 0.0 ) );

            // checking derivative in the 2 direction of nodes 0-3
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 0, 2 ) == Approx( -1.5 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 1, 2 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 2, 2 ) == Approx( 0.0 ) );
            CHECK( tCell.get_cell_info()->compute_cell_size_deriv_general( &tCell, 3, 2 ) == Approx( 1.5 ) );
        }
    }

}    // namespace moris::mtk
