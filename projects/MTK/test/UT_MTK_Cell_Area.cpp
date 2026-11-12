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
#include "fn_norm.hpp"
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
#include "cl_MTK_Cell_ISC.hpp"
#include "cl_MTK_Cell_Proxy.hpp"
#include "cl_MTK_Vertex_Proxy.hpp"
#include "cl_MTK_Cell_Info_Quad9.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Cluster_Proxy.hpp"
#include <limits>
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Interpolation_Function.hpp"

#include <memory>
#include <cmath>
#include "fn_trans.hpp"

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

namespace moris::mtk
{
    namespace
    {
        template< typename Cell_Info_Type >
        class Counted_Cell_Info : public Cell_Info_Type
        {
          public:
            mutable uint mShapeCalls             = 0;
            mutable uint mSpecialCalls           = 0;
            mutable uint mDerivativeCalls        = 0;
            mutable uint mBulkDerivativeCalls    = 0;
            mutable uint mGeneralDerivativeCalls = 0;

            real compute_cell_size_deriv_general( const Cell* aCell, uint aVertex, uint aDirection ) const override
            {
                ++mGeneralDerivativeCalls;

                return Cell_Info_Type::compute_cell_size_deriv_general( aCell, aVertex, aDirection );
            }

            void compute_cell_size_derivatives_general( const Cell* aCell, Matrix< DDRMat >& aDerivatives, const Matrix< DDSMat >& aRequested = {} ) const override
            {
                ++mBulkDerivativeCalls;
                Cell_Info_Type::compute_cell_size_derivatives_general( aCell, aDerivatives, aRequested );
            }

            CellShape compute_cell_shape( const Cell* aCell ) const override
            {
                ++mShapeCalls;

                return Cell_Info_Type::compute_cell_shape( aCell );
            }

            real compute_cell_size_special( const Cell* aCell ) const override
            {
                ++mSpecialCalls;

                return Cell_Info_Type::compute_cell_size_special( aCell );
            }

            real compute_cell_size_deriv( const Cell* aCell, uint aVertexOrdinal, uint aDirection ) const override
            {
                ++mDerivativeCalls;

                return Cell_Info_Type::compute_cell_size_deriv( aCell, aVertexOrdinal, aDirection );
            }
        };

        // Use the normal Cell implementation with owned, mutable test vertices.
        // All vertex addresses remain stable for the lifetime of the cell.
        struct Cell_Measure_Test_Data
        {
            Vector< Vertex_Proxy > mVertices;
            Cell_ISC               mCell;

            static Vector< Vertex* > create_vertices(
                    Vector< Vertex_Proxy >& aVertices,
                    const Matrix< DDRMat >& aCoordinates )
            {
                Vector< Vertex* > tPointers( aVertices.size() );

                for ( uint iVertex = 0; iVertex < aVertices.size(); ++iVertex )
                {
                    aVertices( iVertex ).mVertexCoord = aCoordinates.get_row( iVertex );
                    tPointers( iVertex )              = &aVertices( iVertex );
                }

                return tPointers;
            }

            Cell_Measure_Test_Data( const Matrix< DDRMat >& aCoordinates,
                    const std::shared_ptr< Cell_Info >&     aCellInfo )
                    : mVertices( aCoordinates.n_rows() )
                    , mCell( 0, 0, 0, aCellInfo, create_vertices( mVertices, aCoordinates ) )
            {
            }
        };

        void check_measure_derivatives( Cell_Measure_Test_Data& aData, real aStep, real aTolerance )
        {
            const CellShape tUnperturbedShape   = aData.mCell.get_cell_shape();
            const real      tUnperturbedMeasure = aData.mCell.compute_cell_measure();

            for ( uint iVertex = 0; iVertex < aData.mVertices.size(); ++iVertex )
            {
                for ( uint iDirection = 0; iDirection < aData.mVertices( iVertex ).mVertexCoord.numel(); ++iDirection )
                {
                    INFO( "Vertex=" << iVertex << " direction=" << iDirection << " step=" << aStep );

                    const real tDerivative         = aData.mCell.compute_cell_measure_deriv( iVertex, iDirection );
                    real&      tCoordinate         = aData.mVertices( iVertex ).mVertexCoord( iDirection );
                    const real tOriginalCoordinate = tCoordinate;

                    // Independent vertex perturbations leave the rectangular subspace.
                    // Bypass classification, including its tolerance, for the FD values.
                    tCoordinate = tOriginalCoordinate + aStep;

                    const real tPlus = aData.mCell.compute_cell_measure_general();
                    tCoordinate      = tOriginalCoordinate - aStep;

                    const real tMinus = aData.mCell.compute_cell_measure_general();
                    tCoordinate       = tOriginalCoordinate;

                    CHECK( ( tPlus - tMinus ) / ( 2.0 * aStep ) == Approx( tDerivative ).margin( aTolerance ) );
                }
            }

            // The bypass must not replace the unperturbed cached shape. Coordinates
            // have been restored before resuming normal measure evaluation.
            CHECK( aData.mCell.get_cell_shape() == tUnperturbedShape );
            CHECK( aData.mCell.compute_cell_measure() == Approx( tUnperturbedMeasure ) );
        }

    }    // namespace

    TEST_CASE( "MTK Cell Measure Shape Cache", "[MTK],[Cell_Shape],[Cell_Measure_Cache]" )
    {
        // Two cells share topology information but have different geometry. Repeated
        // shape/measure queries classify each cell once; coordinate and info updates
        // invalidate it. Constant shape does not imply constant area.
        auto                   tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad4 > >();
        Cell_Measure_Test_Data tRectangle( { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } }, tInfo );
        Cell_Measure_Test_Data tSkew( { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.5, 1.0 }, { 0.5, 1.0 } }, tInfo );

        CHECK( tRectangle.mCell.get_cell_shape() == CellShape::RECTANGULAR );
        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tInfo->mShapeCalls == 1 );
        CHECK( tInfo->mSpecialCalls == 2 );
        CHECK( tSkew.mCell.get_cell_shape() == CellShape::PARALLEL );
        CHECK( tSkew.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tInfo->mShapeCalls == 2 );
        CHECK( tInfo->mSpecialCalls == 2 );

        tRectangle.mVertices( 1 ).mVertexCoord( 0 ) = 1.4;
        tRectangle.mCell.invalidate_cell_shape();

        CHECK( tRectangle.mCell.get_cell_shape() == CellShape::STRAIGHT );
        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 1.2 ) );
        CHECK( tInfo->mShapeCalls == 3 );
        CHECK( tInfo->mSpecialCalls == 2 );

        tRectangle.mVertices( 1 ).mVertexCoord( 0 ) = 1.0;
        tRectangle.mCell.invalidate_cell_shape();

        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tInfo->mShapeCalls == 4 );

        auto tReplacementInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad4 > >();
        tRectangle.mCell.set_mtk_cell_info( tReplacementInfo );

        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tReplacementInfo->mShapeCalls == 1 );

        // The proxy uses borrowed Cell_Info rather than the base class shared pointer.
        Cell_Proxy tProxy;
        tProxy.mCellInfo     = tInfo.get();
        tProxy.mVertices     = tRectangle.mCell.get_vertex_pointers();
        tProxy.mSpatialDim   = 2;
        tProxy.mGeometryType = Geometry_Type::QUAD;
        tProxy.mInterpOrder  = Interpolation_Order::LINEAR;
        tProxy.mIntegOrder   = Integration_Order::QUAD_2x2;

        CHECK( tProxy.get_cell_shape() == CellShape::RECTANGULAR );
        CHECK( tProxy.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tInfo->mShapeCalls == 5 );

        // A rectangular geometry update must recompute the measure too; only
        // the classification is cached, never the area.
        tRectangle.mVertices( 1 ).mVertexCoord( 0 ) = 2.0;
        tRectangle.mVertices( 2 ).mVertexCoord( 0 ) = 2.0;
        tRectangle.mCell.invalidate_cell_shape();
        tProxy.invalidate_cell_shape();

        CHECK( tRectangle.mCell.get_cell_shape() == CellShape::RECTANGULAR );
        CHECK( tRectangle.mCell.compute_cell_measure() == Approx( 2.0 ) );
        CHECK( tReplacementInfo->mShapeCalls == 2 );
    }

    TEST_CASE( "MTK Cell Measure Shape Dispatch", "[MTK],[Size_Integration],[Cell_Measure_Cache]" )
    {
        // Rotated rectangles and skew parallelograms have unit area/volume, even
        // though opposite-corner coordinate products give a different result.
        auto                   tQuadInfo        = std::make_shared< Counted_Cell_Info< Cell_Info_Quad4 > >();
        const real             tCos             = std::sqrt( 0.5 );
        Matrix< DDRMat >       tQuadCoordinates = { { 0.0, 0.0 }, { tCos, tCos }, { 0.0, 2.0 * tCos }, { -tCos, tCos } };
        Cell_Measure_Test_Data tRotatedQuad( tQuadCoordinates, tQuadInfo );

        CHECK( tRotatedQuad.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tQuadInfo->mSpecialCalls == 0 );

        auto             tHexInfo        = std::make_shared< Counted_Cell_Info< Cell_Info_Hex8 > >();
        Matrix< DDRMat > tHexCoordinates = {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }
        };

        Cell_Measure_Test_Data tBox( tHexCoordinates, tHexInfo );

        CHECK( tBox.mCell.compute_cell_measure_general() == Approx( 1.0 ) );
        CHECK( tHexInfo->mShapeCalls == 0 );
        CHECK( tHexInfo->mSpecialCalls == 0 );
        CHECK( tBox.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tHexInfo->mSpecialCalls == 1 );
        CHECK( tHexInfo->mShapeCalls == 1 );

        Matrix< DDRMat > tRotatedCoordinates = tHexCoordinates;

        for ( uint iVertex = 0; iVertex < 8; ++iVertex )
        {
            tRotatedCoordinates( iVertex, 0 ) = tCos * ( tHexCoordinates( iVertex, 0 ) - tHexCoordinates( iVertex, 1 ) );
            tRotatedCoordinates( iVertex, 1 ) = tCos * ( tHexCoordinates( iVertex, 0 ) + tHexCoordinates( iVertex, 1 ) );
        }

        Cell_Measure_Test_Data tRotatedHex( tRotatedCoordinates, tHexInfo );

        CHECK( tRotatedHex.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tHexInfo->mSpecialCalls == 1 );

        for ( uint iVertex = 4; iVertex < 8; ++iVertex )
        {
            tHexCoordinates( iVertex, 0 ) += 0.5;
        }

        Cell_Measure_Test_Data tSkewHex( tHexCoordinates, tHexInfo );

        CHECK( tSkewHex.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tHexInfo->mSpecialCalls == 1 );

        // A single lifted corner gives trilinear volume 1 + 0.4/4.
        tBox.mVertices( 6 ).mVertexCoord( 2 ) = 1.4;
        tBox.mCell.invalidate_cell_shape();

        CHECK( tBox.mCell.compute_cell_measure() == Approx( 1.1 ) );
        CHECK( tHexInfo->mSpecialCalls == 1 );
    }

    TEST_CASE( "MTK Cached Quad4 Measure Sensitivity", "[MTK],[Size_Deriv],[Cell_Measure_Cache]" )
    {
        // Cache a unit square, then perturb every vertex coordinate independently.
        // Moving vertex 1 in x has dA/dx=1/2; the rectangular shortcut would give zero.
        auto                   tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad4 > >();
        Cell_Measure_Test_Data tData( { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } }, tInfo );

        CHECK( tData.mCell.compute_cell_measure() == Approx( 1.0 ) );
        CHECK( tData.mCell.compute_cell_measure_deriv( 1, 0 ) == Approx( 0.5 ) );
        check_measure_derivatives( tData, 1e-5, 1e-8 );
        check_measure_derivatives( tData, 1e-10, 5e-6 );

        CHECK( tInfo->mShapeCalls == 1 );
    }

    TEST_CASE( "MTK Cached Hex8 Measure Sensitivity", "[MTK],[Size_Deriv],[Cell_Measure_Cache],[Hex8_Measure_Deriv]" )
    {
        // Each coordinate of a unit cube corner has volume derivative +/-1/4.
        // FD steps both above and below the shape tolerance use general integration.
        auto             tInfo        = std::make_shared< Counted_Cell_Info< Cell_Info_Hex8 > >();
        Matrix< DDRMat > tCoordinates = {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }
        };

        Cell_Measure_Test_Data tData( tCoordinates, tInfo );

        CHECK( tData.mCell.compute_cell_measure() == Approx( 1.0 ) );

        for ( uint iVertex = 0; iVertex < 8; ++iVertex )
        {
            for ( uint iDirection = 0; iDirection < 3; ++iDirection )
            {
                CHECK( tData.mCell.compute_cell_measure_deriv( iVertex, iDirection )
                        == Approx( 0.25 * ( 2.0 * tCoordinates( iVertex, iDirection ) - 1.0 ) ) );
            }
        }

        check_measure_derivatives( tData, 1e-5, 1e-8 );
        check_measure_derivatives( tData, 1e-10, 5e-6 );

        CHECK( tInfo->mShapeCalls == 1 );
    }

    TEST_CASE( "MTK Analytical Hex8 Measure Derivatives", "[MTK],[Size_Deriv],[Hex8_Measure_Deriv]" )
    {
        // Compare every analytical nodal derivative with integration and FD on
        // boxes, rotated/skewed affine cells, and a non-affine cell with warped
        // faces. Translation invariance and the cubic volume scaling provide
        // independent checks of all 24 derivatives.
        const Matrix< DDRMat > tUnitCoordinates = {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }
        };

        auto tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Hex8 > >();

        for ( uint iCase = 0; iCase < 4; ++iCase )
        {
            DYNAMIC_SECTION( "HEX8 geometry=" << iCase )
            {
                Matrix< DDRMat > tCoordinates = tUnitCoordinates;

                for ( uint iVertex = 0; iVertex < 8; ++iVertex )
                {
                    const real tX = tUnitCoordinates( iVertex, 0 );
                    const real tY = tUnitCoordinates( iVertex, 1 );
                    const real tZ = tUnitCoordinates( iVertex, 2 );

                    if ( iCase == 1 )
                    {
                        tCoordinates( iVertex, 0 ) = 2.0 * tX;
                        tCoordinates( iVertex, 1 ) = 3.0 * tY;
                        tCoordinates( iVertex, 2 ) = 4.0 * tZ;
                    }
                    else if ( iCase == 2 )
                    {
                        const real tCos            = std::sqrt( 0.5 );
                        tCoordinates( iVertex, 0 ) = 17.0 + tCos * ( 2.0 * tX - 3.0 * tY ) + 0.4 * tZ;
                        tCoordinates( iVertex, 1 ) = -23.0 + tCos * ( 2.0 * tX + 3.0 * tY ) + 0.3 * tZ;
                        tCoordinates( iVertex, 2 ) = 4.0 + 4.0 * tZ;
                    }
                }

                if ( iCase == 3 )
                {
                    tCoordinates( 6, 2 ) += 0.7;
                    tCoordinates( 1, 1 ) -= 0.2;
                    tCoordinates( 3, 0 ) -= 0.15;
                    tCoordinates( 7, 1 ) += 0.25;
                }

                Cell_Measure_Test_Data tData( tCoordinates, tInfo );
                Matrix< DDRMat >       tDerivativeSum( 3, 1, 0.0 );
                real                   tScalingDerivative = 0.0;
                const uint             tInitialCalls      = tInfo->mDerivativeCalls;

                for ( uint iVertex = 0; iVertex < 8; ++iVertex )
                {
                    for ( uint iDirection = 0; iDirection < 3; ++iDirection )
                    {
                        INFO( "Vertex=" << iVertex << " direction=" << iDirection );

                        const real tDerivative = tData.mCell.compute_cell_measure_deriv( iVertex, iDirection );

                        CHECK( tDerivative == Approx( tInfo->compute_cell_size_deriv_general( &tData.mCell, iVertex, iDirection ) ).epsilon( 1e-11 ).margin( 1e-10 ) );
                        tDerivativeSum( iDirection ) += tDerivative;
                        tScalingDerivative += tCoordinates( iVertex, iDirection ) * tDerivative;
                    }
                }

                CHECK( tInfo->mDerivativeCalls - tInitialCalls == 24 );

                for ( uint iDirection = 0; iDirection < 3; ++iDirection )
                {
                    CHECK( tDerivativeSum( iDirection ) == Approx( 0.0 ).margin( 1e-12 ) );
                }

                CHECK( tScalingDerivative == Approx( 3.0 * tData.mCell.compute_cell_measure_general() ).epsilon( 1e-11 ).margin( 1e-10 ) );
                check_measure_derivatives( tData, 1e-5, 1e-7 );
            }
        }
    }

    TEST_CASE( "MTK Curved Quad9 Measure Sensitivity", "[MTK],[Size_Deriv],[Cell_Measure_Cache]" )
    {
        // Bow the top edge of [-1,1]^2 with y=eta+0.3*(1-xi^2)*(1+eta)/2.
        // Integrating the boundary gives area 4.4. Midside sensitivities must use
        // the full isoparametric geometry, not a formula using only corner nodes.
        auto                                           tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad9 > >();
        Interpolation_Rule                             tRule( Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, Interpolation_Order::QUADRATIC, Interpolation_Type::LAGRANGE, Interpolation_Order::LINEAR );
        std::unique_ptr< Interpolation_Function_Base > tFunction( tRule.create_space_interpolation_function() );

        Matrix< DDRMat > tParametricCoordinates;
        tFunction->get_param_coords( tParametricCoordinates );

        Matrix< DDRMat >       tCoordinates = trans( tParametricCoordinates );
        Cell_Measure_Test_Data tData( tCoordinates, tInfo );

        CHECK( tData.mCell.compute_cell_measure() == Approx( 4.0 ) );

        for ( uint iVertex = 0; iVertex < tCoordinates.n_rows(); ++iVertex )
        {
            const real tXi  = tCoordinates( iVertex, 0 );
            const real tEta = tCoordinates( iVertex, 1 );
            tData.mVertices( iVertex ).mVertexCoord( 1 ) += 0.3 * ( 1.0 - tXi * tXi ) * ( 1.0 + tEta ) / 2.0;
        }

        tData.mCell.invalidate_cell_shape();

        CHECK( tData.mCell.get_cell_shape() == CellShape::GENERAL );
        CHECK( tData.mCell.compute_cell_measure() == Approx( 4.4 ) );
        CHECK( tInfo->mSpecialCalls == 1 );
        check_measure_derivatives( tData, 1e-5, 1e-8 );

        CHECK( tInfo->mShapeCalls == 2 );
    }

    TEST_CASE( "MTK Linear Simplex Measure Shortcuts", "[MTK],[Size_Integration],[Cell_Measure_Cache]" )
    {
        // Arbitrary nonrectangular linear simplices retain their determinant-based
        // shortcuts. Unit triangle and tetrahedron references are 1/2 and 1/6.
        auto                   tTriInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Tri3 > >();
        Cell_Measure_Test_Data tTriangle( { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } }, tTriInfo );

        CHECK( tTriangle.mCell.compute_cell_measure() == Approx( 0.5 ) );
        CHECK( tTriInfo->mSpecialCalls == 1 );
        CHECK( tTriInfo->mShapeCalls == 0 );
        check_measure_derivatives( tTriangle, 1e-5, 1e-8 );

        auto                   tTetInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Tet4 > >();
        Cell_Measure_Test_Data tTetrahedron( { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } }, tTetInfo );

        CHECK( tTetrahedron.mCell.compute_cell_measure() == Approx( 1.0 / 6.0 ) );
        CHECK( tTetInfo->mSpecialCalls == 1 );
        CHECK( tTetInfo->mShapeCalls == 0 );
        check_measure_derivatives( tTetrahedron, 1e-5, 1e-8 );
    }

    TEST_CASE( "MTK Cell Shape Scale Invariance", "[MTK],[Cell_Shape],[Cell_Shape_Scale]" )
    {
        // Identical mappings must retain their classification and normalized
        // measure across scales. In particular, the tiny trapezoid previously
        // selected the rectangular shortcut and overestimated its area by 20%.
        Cell_Info_Factory tFactory;

        for ( CellTopology tTopology : { CellTopology::QUAD4, CellTopology::QUAD8, CellTopology::QUAD9, CellTopology::QUAD16, CellTopology::HEX8, CellTopology::HEX20, CellTopology::HEX27, CellTopology::HEX64 } )
        {
            auto             tInfo = tFactory.create_cell_info_sp( tTopology );
            Matrix< DDRMat > tReference;
            tInfo->get_loc_coords_of_cell( tReference );

            const uint tDimension   = tReference.n_cols();
            const bool tHigherOrder = tInfo->get_num_verts() > ( tDimension == 2 ? 4u : 8u );

            for ( real tScale : { 1.0e-9, 1.0e-5, 1.0, 1.0e6 } )
            {
                for ( bool tTranslate : { false, true } )
                {
                    for ( uint iMapping = 0; iMapping < ( tHigherOrder ? 6u : 5u ); ++iMapping )
                    {
                        if ( tDimension == 2 && iMapping == 4 )
                        {
                            continue;
                        }

                        DYNAMIC_SECTION( "Topology=" << static_cast< uint >( tTopology ) << " scale=" << tScale
                                                     << " translated=" << tTranslate << " mapping=" << iMapping )
                        {
                            Matrix< DDRMat > tCoordinates( tReference.n_rows(), tDimension );

                            for ( uint iNode = 0; iNode < tCoordinates.n_rows(); ++iNode )
                            {
                                const real tXi = tReference( iNode, 0 );
                                const real tU  = ( tXi + 1.0 ) / 2.0;
                                const real tV  = ( tReference( iNode, 1 ) + 1.0 ) / 2.0;
                                const real tW  = tDimension == 3 ? ( tReference( iNode, 2 ) + 1.0 ) / 2.0 : 0.0;
                                real       tX  = tU;
                                real       tY  = tV;
                                real       tZ  = tW;

                                if ( iMapping == 1 )
                                {
                                    tX = 0.8 * tU - 0.6 * tV;
                                    tY = 0.6 * tU + 0.8 * tV;
                                }
                                else if ( iMapping == 2 )
                                {
                                    tX += 0.3 * tV;
                                }
                                else if ( iMapping == 3 )
                                {
                                    // A trapezoid or a frustum with planar faces.
                                    tX *= 1.0 + 0.5 * ( tDimension == 2 ? tV : tW );

                                    if ( tDimension == 3 )
                                    {
                                        tY *= 1.0 + 0.5 * tW;
                                    }
                                }
                                else if ( iMapping == 4 )
                                {
                                    tZ += 0.2 * tU * tV * tW;
                                }
                                else if ( iMapping == 5 )
                                {
                                    if ( tDimension == 2 )
                                    {
                                        tY += 0.3 * ( 1.0 - tXi * tXi ) * tV;
                                    }
                                    else
                                    {
                                        tZ += 0.3 * ( 1.0 - tXi * tXi ) * tV * tW;
                                    }
                                }

                                tCoordinates( iNode, 0 ) = tScale * ( tX + ( tTranslate ? 8.0 : 0.0 ) );
                                tCoordinates( iNode, 1 ) = tScale * ( tY - ( tTranslate ? 4.0 : 0.0 ) );

                                if ( tDimension == 3 )
                                {
                                    tCoordinates( iNode, 2 ) = tScale * ( tZ + ( tTranslate ? 2.0 : 0.0 ) );
                                }
                            }

                            Cell_Measure_Test_Data tData( tCoordinates, tInfo );
                            const CellShape        tExpected = iMapping == 0 ? CellShape::RECTANGULAR
                                                             : iMapping < 3  ? CellShape::PARALLEL
                                                             : iMapping == 3 ? CellShape::STRAIGHT
                                                                             : CellShape::GENERAL;

                            CHECK( tData.mCell.get_cell_shape() == tExpected );

                            const real tMeasureScale = std::pow( tScale, tDimension );
                            const real tMeasure      = tData.mCell.compute_cell_measure() / tMeasureScale;
                            const real tIntegrated   = tData.mCell.compute_cell_measure_general() / tMeasureScale;

                            CHECK( tMeasure == Approx( tIntegrated ).epsilon( 1e-10 ) );

                            if ( iMapping == 3 )
                            {
                                CHECK( tMeasure == Approx( tDimension == 2 ? 1.25 : 19.0 / 12.0 ).epsilon( 1e-10 ) );
                            }
                            else if ( iMapping < 3 )
                            {
                                CHECK( tMeasure == Approx( 1.0 ).epsilon( 1e-10 ) );
                            }
                        }
                    }
                }
            }
        }
    }

    TEST_CASE( "MTK Cell Shape Degenerate Geometry", "[MTK],[Cell_Shape],[Cell_Shape_Scale]" )
    {
        // Degenerate or nonfinite directions must never select the rectangular
        // shortcut. Exercise every classifier without integrating invalid cells.
        Cell_Info_Factory tFactory;

        for ( CellTopology tTopology : { CellTopology::QUAD4, CellTopology::QUAD8, CellTopology::QUAD9, CellTopology::QUAD16, CellTopology::HEX8, CellTopology::HEX20, CellTopology::HEX27, CellTopology::HEX64 } )
        {
            auto             tInfo = tFactory.create_cell_info_sp( tTopology );
            Matrix< DDRMat > tReference;
            tInfo->get_loc_coords_of_cell( tReference );

            for ( uint iInvalid = 0; iInvalid < 3; ++iInvalid )
            {
                Matrix< DDRMat > tCoordinates = tReference;

                for ( uint iNode = 0; iNode < tCoordinates.n_rows(); ++iNode )
                {
                    tCoordinates( iNode, 1 ) = 0.0;
                }

                if ( iInvalid > 0 )
                {
                    tCoordinates         = tReference;
                    tCoordinates( 0, 0 ) = iInvalid == 1 ? std::numeric_limits< real >::quiet_NaN() : std::numeric_limits< real >::infinity();
                }

                Cell_Measure_Test_Data tData( tCoordinates, tInfo );

                CHECK( tData.mCell.get_cell_shape() == CellShape::GENERAL );
            }

            // Extreme aspect ratios retain their axis-aligned classification;
            // each direction must be normalized using its own length.
            for ( real tAspect : { 1.0e-12, 1.0e12 } )
            {
                Matrix< DDRMat > tCoordinates = tReference;

                for ( uint iNode = 0; iNode < tCoordinates.n_rows(); ++iNode )
                {
                    tCoordinates( iNode, 0 ) *= tAspect;
                }

                Cell_Measure_Test_Data tData( tCoordinates, tInfo );

                CHECK( tData.mCell.get_cell_shape() == CellShape::RECTANGULAR );
            }
        }
    }

    TEST_CASE( "MTK Tiny Cell Shape Cache Invalidation", "[MTK],[Cell_Shape_Scale],[Cell_Measure_Cache]" )
    {
        // Classification stays lazy, but an invalidated tiny square must switch
        // to general integration after its upper-right corner becomes skewed.
        auto                   tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad4 > >();
        Cell_Measure_Test_Data tData( { { 0.0, 0.0 }, { 1e-5, 0.0 }, { 1e-5, 1e-5 }, { 0.0, 1e-5 } }, tInfo );

        CHECK( tData.mCell.compute_cell_measure() / 1e-10 == Approx( 1.0 ) );
        CHECK( tData.mCell.get_cell_shape() == CellShape::RECTANGULAR );
        CHECK( tInfo->mShapeCalls == 1 );
        tData.mVertices( 2 ).mVertexCoord( 0 ) = 1.5e-5;
        tData.mCell.invalidate_cell_shape();

        CHECK( tData.mCell.get_cell_shape() == CellShape::STRAIGHT );
        CHECK( tData.mCell.compute_cell_measure() / 1e-10 == Approx( 1.25 ) );
        CHECK( tInfo->mShapeCalls == 2 );
        CHECK( tInfo->mSpecialCalls == 1 );
    }

    TEST_CASE( "MTK Bulk Cell Measure Derivatives", "[MTK],[Bulk_Measure_Derivatives]" )
    {
        // Warped polynomial maps exercise all supported linear and higher-order
        // cells, including simplex barycentric normalization and midside nodes.
        Cell_Info_Factory tFactory;

        for ( CellTopology tTopology : { CellTopology::TRI3, CellTopology::TRI6, CellTopology::TET4, CellTopology::TET10, CellTopology::QUAD4, CellTopology::QUAD8, CellTopology::QUAD9, CellTopology::QUAD16, CellTopology::HEX8, CellTopology::HEX20, CellTopology::HEX27, CellTopology::HEX64 } )
        {
            INFO( "Topology=" << static_cast< uint >( tTopology ) );

            auto             tInfo = tFactory.create_cell_info_sp( tTopology );
            Matrix< DDRMat > tReference;
            tInfo->get_loc_coords_of_cell( tReference );

            Matrix< DDRMat > tCoordinates = tReference;

            for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
            {
                const real tXi       = tReference( i, 0 );
                const real tEta      = tReference( i, 1 );
                tCoordinates( i, 0 ) = 1.2 * tXi + 0.2 * tEta + 0.05 * tXi * tEta;
                tCoordinates( i, 1 ) = 0.1 * tXi + 1.1 * tEta + 0.08 * tXi * tXi;

                if ( tCoordinates.n_cols() == 3 )
                {
                    const real tZeta     = tReference( i, 2 );
                    tCoordinates( i, 2 ) = 0.1 * tXi + 0.2 * tEta + 1.3 * tZeta + 0.06 * tXi * tEta + 0.04 * tXi * tZeta;
                }
            }

            Cell_Measure_Test_Data tData( tCoordinates, tInfo );
            Matrix< DDRMat >       tBulk;
            tData.mCell.compute_cell_measure_derivatives( tBulk );

            REQUIRE( tBulk.n_rows() == tCoordinates.n_rows() );
            REQUIRE( tBulk.n_cols() == tCoordinates.n_cols() );

            for ( uint i = 0; i < tBulk.n_rows(); ++i )
            {
                for ( uint d = 0; d < tBulk.n_cols(); ++d )
                {
                    INFO( "Vertex=" << i << " direction=" << d );
                    CHECK( tBulk( i, d ) == Approx( tData.mCell.compute_cell_measure_deriv( i, d ) ).margin( 1e-11 ) );
                    CHECK( tBulk( i, d ) == Approx( tInfo->compute_cell_size_deriv_general( &tData.mCell, i, d ) ).margin( 1e-11 ) );

                    real&      tCoordinate = tData.mVertices( i ).mVertexCoord( d );
                    const real tOriginal   = tCoordinate;

                    // Subtracting two integrated volumes amplifies quadrature
                    // roundoff at tiny steps, especially for HEX64's zero
                    // derivatives. Check two larger steps without relaxing the
                    // tolerance; these maps retain their Jacobian orientation.
                    for ( real tStep : { 1e-4, 5e-5 } )
                    {
                        INFO( "Finite-difference step=" << tStep );
                        tCoordinate = tOriginal + tStep;

                        const real tPlus = tData.mCell.compute_cell_measure_general();
                        tCoordinate      = tOriginal - tStep;

                        const real tMinus = tData.mCell.compute_cell_measure_general();
                        tCoordinate       = tOriginal;

                        CHECK( tBulk( i, d ) == Approx( ( tPlus - tMinus ) / ( 2.0 * tStep ) ).margin( 1e-8 ) );
                    }
                }
            }

            // Restrict directions and vertices, then change geometry and repeat:
            // neither the scalar derivative cache nor a geometry cache may leak in.
            Matrix< DDSMat > tMask( tBulk.n_rows(), tBulk.n_cols(), 0 );
            tMask( 0, 0 )                                   = 1;
            tMask( tBulk.n_rows() - 1, tBulk.n_cols() - 1 ) = 1;

            Matrix< DDRMat > tPartial;
            tData.mCell.compute_cell_measure_derivatives( tPartial, tMask );

            for ( uint i = 0; i < tBulk.n_rows(); ++i )
            {
                for ( uint d = 0; d < tBulk.n_cols(); ++d )
                {
                    CHECK( tPartial( i, d ) == Approx( tMask( i, d ) ? tBulk( i, d ) : 0.0 ).margin( 1e-12 ) );
                }
            }

            tData.mVertices( 1 ).mVertexCoord( 1 ) += 0.03;
            tData.mCell.compute_cell_measure_derivatives( tPartial, tMask );

            CHECK( tPartial( 0, 0 ) == Approx( tData.mCell.compute_cell_measure_deriv( 0, 0 ) ).margin( 1e-11 ) );
            tData.mVertices( 1 ).mVertexCoord( 1 ) = tCoordinates( 1, 1 );
            tMask.fill( 0 );
            tData.mCell.compute_cell_measure_derivatives( tPartial, tMask );

            CHECK( norm( tPartial ) == 0.0 );

            Matrix< DDSMat > tWrongMask( 1, 1, 1 );

            REQUIRE_THROWS( tData.mCell.compute_cell_measure_derivatives( tPartial, tWrongMask ) );
        }
    }

    TEST_CASE( "MTK Bulk Cluster Measure Derivatives", "[MTK],[Bulk_Measure_Derivatives]" )
    {
        // Two curved QUAD9 cells share an edge. Duplicate requested coordinates
        // preserve the scalar API's enriched-node correspondence; unmatched and
        // inactive coordinates must contribute zero for primary and void cells.
        auto             tInfo = std::make_shared< Counted_Cell_Info< Cell_Info_Quad9 > >();
        Matrix< DDRMat > tCoordinates;
        tInfo->get_loc_coords_of_cell( tCoordinates );

        for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
        {
            tCoordinates( i, 1 ) += 0.2 * ( 1.0 - tCoordinates( i, 0 ) * tCoordinates( i, 0 ) );
        }

        Cell_Measure_Test_Data tFirst( tCoordinates, tInfo );
        Matrix< DDRMat >       tSecondCoordinates = tCoordinates;

        for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
        {
            tSecondCoordinates( i, 0 ) += 2.0;
        }

        Cell_Measure_Test_Data    tSecond( tSecondCoordinates, tInfo );
        moris::Cell_Cluster_Proxy tCluster;
        tCluster.mPrimaryIntegrationCells = { &tFirst.mCell, &tSecond.mCell };
        tCluster.mVoidIntegrationCells    = { &tSecond.mCell };

        Matrix< DDRMat > tQueries( 20, 2, 0.0 );

        for ( uint i = 0; i < 9; ++i )
        {
            tQueries.get_row( i )     = tCoordinates.get_row( i );
            tQueries.get_row( i + 9 ) = tSecondCoordinates.get_row( i );
        }

        tQueries.get_row( 18 ) = tCoordinates.get_row( 0 );
        tQueries( 19, 0 )      = 99.0;

        Matrix< DDSMat > tMask( 20, 2, 1 );
        tMask( 4, 0 ) = 0;

        for ( Primary_Void tPhase : { Primary_Void::PRIMARY, Primary_Void::VOID } )
        {
            Matrix< DDRMat > tExpected( 20, 2, 0.0 );

            for ( uint i = 0; i < 20; ++i )
            {
                for ( uint d = 0; d < 2; ++d )
                {
                    if ( tMask( i, d ) )
                    {
                        tExpected( i, d ) = tCluster.compute_cluster_cell_measure_derivative( tQueries.get_row( i ), d, tPhase );
                    }
                }
            }

            tInfo->mBulkDerivativeCalls    = 0;
            tInfo->mGeneralDerivativeCalls = 0;

            const Matrix< DDRMat > tActual = tCluster.compute_cluster_cell_measure_derivatives( tQueries, tMask, tPhase );

            CHECK( norm( tActual - tExpected ) < 1e-11 );
            CHECK( tInfo->mBulkDerivativeCalls == ( tPhase == Primary_Void::PRIMARY ? 2 : 1 ) );
            CHECK( tInfo->mGeneralDerivativeCalls == 0 );
        }

        tInfo->mBulkDerivativeCalls    = 0;
        tInfo->mGeneralDerivativeCalls = 0;
        tMask.fill( 0 );

        Matrix< DDRMat > tActual = tCluster.compute_cluster_cell_measure_derivatives( tQueries, tMask );

        CHECK( norm( tActual ) == 0.0 );
        CHECK( tInfo->mBulkDerivativeCalls == 0 );
        tMask( 0, 0 ) = 1;
        tActual       = tCluster.compute_cluster_cell_measure_derivatives( tQueries, tMask );

        CHECK( tInfo->mBulkDerivativeCalls == 0 );
        CHECK( tInfo->mGeneralDerivativeCalls == 1 );
        CHECK( tActual( 0, 0 ) == Approx( tFirst.mCell.compute_cell_measure_deriv( 0, 0 ) ) );
        tMask.fill( 0 );
        tMask( 19, 0 ) = 1;
        tActual        = tCluster.compute_cluster_cell_measure_derivatives( tQueries, tMask );

        CHECK( norm( tActual ) == 0.0 );
    }

    TEST_CASE( "MTK Bulk Determinant Derivative Cache Independence", "[MTK],[Bulk_Measure_Derivatives]" )
    {
        // Reuse one interpolator at multiple points, with its scalar derivative
        // cache already populated. Bulk evaluation must not reuse that scalar.
        auto             tInfo = std::make_shared< Cell_Info_Quad9 >();
        Matrix< DDRMat > tCoordinates;
        tInfo->get_loc_coords_of_cell( tCoordinates );
        tCoordinates( 6, 1 ) += 0.3;

        Interpolation_Rule tRule( Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, Interpolation_Order::QUADRATIC, Interpolation_Type::LAGRANGE, Interpolation_Order::LINEAR );
        Space_Interpolator tInterpolator( tRule );
        tInterpolator.set_space_coeff( tCoordinates );
        tInterpolator.set_param_coeff();

        for ( real tXi : { -0.4, 0.7 } )
        {
            tInterpolator.set_space( Matrix< DDRMat >{ { tXi }, { 0.2 } } );

            const real       tFirst = tInterpolator.space_det_J_deriv( 0, 0 );
            Matrix< DDRMat > tBulk;
            tInterpolator.space_det_J_derivatives( tBulk );

            CHECK( tInterpolator.space_det_J_deriv( 0, 0 ) == tFirst );

            for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
            {
                for ( uint d = 0; d < 2; ++d )
                {
                    tInterpolator.reset_eval_flags_deriv();

                    CHECK( tBulk( i, d ) == Approx( tInterpolator.space_det_J_deriv( i, d ) ).margin( 1e-12 ) );
                }
            }
        }
    }

    TEST_CASE( "MTK Bulk Cluster Scalar Fallback", "[MTK],[Bulk_Measure_Derivatives]" )
    {
        // Cluster implementations without a bulk specialization must retain their
        // scalar virtual behavior and only receive requested directions.
        class Scalar_Cluster : public moris::Cell_Cluster_Proxy
        {
          public:
            mutable uint mCalls = 0;
            real         compute_cluster_cell_measure_derivative( const Matrix< DDRMat >& aCoordinates, uint aDirection, Primary_Void aPhase, Leader_Follower aSide ) const override
            {
                ++mCalls;

                return aCoordinates( 0 ) + aDirection + ( aPhase == Primary_Void::VOID ? 10.0 : 0.0 )
                     + ( aSide == Leader_Follower::FOLLOWER ? 20.0 : 0.0 );
            }
        };

        Scalar_Cluster         tCluster;
        Matrix< DDRMat >       tCoordinates = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        Matrix< DDSMat >       tRequested   = { { 1, 0 }, { 0, 1 } };
        const Matrix< DDRMat > tResult      = tCluster.Cluster::compute_cluster_cell_measure_derivatives(
                tCoordinates, tRequested, Primary_Void::VOID, Leader_Follower::FOLLOWER );

        CHECK( tCluster.mCalls == 2 );
        CHECK( tResult( 0, 0 ) == Approx( 31.0 ) );
        CHECK( tResult( 1, 1 ) == Approx( 34.0 ) );
        CHECK( tResult( 0, 1 ) == 0.0 );
        CHECK( tResult( 1, 0 ) == 0.0 );
    }

}    // namespace moris::mtk
