/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Hex20.cpp
 *
 */

#include "cl_MTK_Cell_Info_Hex20.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------
        enum Geometry_Type
        Cell_Info_Hex20::get_cell_geometry() const
        {
            return Geometry_Type::HEX;
        }

        // ----------------------------------------------------------------------------------
        enum CellTopology
        Cell_Info_Hex20::get_cell_topology() const
        {
            return CellTopology::HEX20;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Hex20::get_cell_interpolation_order() const
        {
            return Interpolation_Order::SERENDIPITY;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Hex20::get_cell_integration_order() const
        {
            return Integration_Order::HEX_3x3x3;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Hex20::compute_cell_shape( moris::mtk::Cell const * aCell ) const
        {
            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Vector< Vertex* > tVertices = aCell->get_vertex_pointers();

            // init cell shape
            CellShape tCellShape = CellShape::RECTANGULAR;

            // error threshold
            real tEpsilon = 1.0E-8;

            // init cell of face normals
            moris::Vector< moris::Matrix< DDRMat > > tFaceNormals( 6 );

            // looping through each face
            for ( uint iFace = 0; iFace < 6; iFace++ )
            {
                // getting nodes on the face
                moris::Matrix< moris::IndexMat > tFaceNodes = this->get_node_to_face_map( iFace );

                // get the vertex coords
                moris::Matrix< DDRMat > tVertex0 = tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex1 = tVertices( tFaceNodes( 1 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex2 = tVertices( tFaceNodes( 2 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex3 = tVertices( tFaceNodes( 3 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex4 = tVertices( tFaceNodes( 4 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex5 = tVertices( tFaceNodes( 5 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex6 = tVertices( tFaceNodes( 6 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex7 = tVertices( tFaceNodes( 7 ) )->get_coords();

                // get edges to define check plane
                auto tEdge0           = tVertex1 - tVertex0;
                auto tEdge1           = tVertex2 - tVertex1;
                tFaceNormals( iFace ) = cross( tEdge0, tEdge1 );

                // these are the other plane normals produced by the other 6 nodes.
                // keeping node 0 to ensure the plane is at the same offset.
                auto tEdge2     = tVertex3 - tVertex0;
                auto tEdge3     = tVertex4 - tVertex3;
                auto tFaceVec23 = cross( tEdge2, tEdge3 );

                auto tEdge4     = tVertex5 - tVertex0;
                auto tEdge5     = tVertex6 - tVertex5;
                auto tFaceVec45 = cross( tEdge4, tEdge5 );

                auto tEdge6     = tVertex7 - tVertex0;
                auto tEdge7     = tVertex4 - tVertex7;
                auto tFaceVec67 = cross( tEdge6, tEdge7 );

                // All three of the plane normals must be parallel in order to be considered straight shape
                if (                                                                        //
                        norm( cross( tFaceNormals( iFace ), tFaceVec23 ) ) > tEpsilon ||    //
                        norm( cross( tFaceNormals( iFace ), tFaceVec45 ) ) > tEpsilon ||    //
                        norm( cross( tFaceNormals( iFace ), tFaceVec67 ) ) > tEpsilon )
                {
                    tCellShape = CellShape::GENERAL;
                    break;
                }
            }

            // since the faces are planar, checking if it is rectangular (aligned), parallel, or straight

            // if the cell shape is straight, determine if it is a parallel cell
            if ( tCellShape == CellShape::RECTANGULAR )
            {
                // create the cross products of opposite faces
                auto tCross02 = cross( tFaceNormals( 0 ), tFaceNormals( 2 ) );
                auto tCross13 = cross( tFaceNormals( 1 ), tFaceNormals( 3 ) );
                auto tCross45 = cross( tFaceNormals( 4 ), tFaceNormals( 5 ) );

                // checking if these face normals are parallel
                if (                                      //
                        norm( tCross02 ) > tEpsilon ||    //
                        norm( tCross13 ) > tEpsilon ||    //
                        norm( tCross45 ) > tEpsilon )
                {
                    tCellShape = CellShape::STRAIGHT;
                }
                else
                {
                    // get the some edges
                    moris::Matrix< DDRMat > tVertex1 = tVertices( 1 )->get_coords();

                    moris::Matrix< DDRMat > tEdge0 = tVertex1 - tVertices( 0 )->get_coords();
                    moris::Matrix< DDRMat > tEdge1 = tVertices( 2 )->get_coords() - tVertex1;

                    // checking that this cell fits requirements for being rectangular and aligned (CellShape::RECTANGULAR)
                    if (                                                                                   //
                            std::abs( dot( tFaceNormals( 0 ), tFaceNormals( 1 ) ) ) > tEpsilon ||          //
                            std::abs( dot( tFaceNormals( 0 ), tFaceNormals( 4 ) ) ) > tEpsilon ||          //
                            std::abs( dot( tFaceNormals( 1 ), tFaceNormals( 4 ) ) ) > tEpsilon ||          //
                            std::abs( tEdge0( 1 ) ) > tEpsilon || std::abs( tEdge0( 2 ) ) > tEpsilon ||    //
                            std::abs( tEdge1( 0 ) ) > tEpsilon || std::abs( tEdge1( 2 ) ) > tEpsilon ||    //
                            tEdge0( 0 ) < 0.0 || tEdge1( 1 ) < 0.0 )
                    {
                        tCellShape = CellShape::PARALLEL;
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex20::get_num_verts() const
        {
            return 20;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex20::get_num_facets() const
        {
            return 6;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex20::get_num_edges() const
        {
            return 12;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex20::get_num_verts_per_facet() const
        {
            return 8;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex20::get_loc_coord_dim() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        inline moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_face_map() const
        {
            return { { 0, 1, 5, 4, 8, 13, 16, 12 },
                { 1, 2, 6, 5, 9, 14, 17, 13 },
                { 2, 3, 7, 6, 10, 15, 18, 14 },
                { 0, 4, 7, 3, 12, 19, 15, 11 },
                { 0, 3, 2, 1, 11, 10, 9, 8 },
                { 4, 5, 6, 7, 16, 17, 18, 19 } };
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_edge_map() const
        {
            return { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        }

        // ----------------------------------------------------------------------------------

        inline moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_face_map( moris::uint aSideOrdinal ) const
        {
            switch ( aSideOrdinal )
            {
                case ( 0 ):
                {
                    return { { 0, 1, 5, 4, 8, 13, 16, 12 } };
                    break;
                }
                case ( 1 ):
                {
                    return { { 1, 2, 6, 5, 9, 14, 17, 13 } };
                    break;
                }
                case ( 2 ):
                {
                    return { { 2, 3, 7, 6, 10, 15, 18, 14 } };
                    break;
                }
                case ( 3 ):
                {
                    return { { 0, 4, 7, 3, 12, 19, 15, 11 } };
                    break;
                }
                case ( 4 ):
                {
                    return { { 0, 3, 2, 1, 11, 10, 9, 8 } };
                    break;
                }
                case ( 5 ):
                {
                    return { { 4, 5, 6, 7, 16, 17, 18, 19 } };
                    break;
                }
                default:
                    MORIS_ERROR( 0, "Invalid side ordinal specified" );
                    return moris::Matrix< moris::IndexMat >( 0, 0 );
                    break;
            }
        }

        // ----------------------------------------------------------------------------------
        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_edge_map( moris::uint aEdgeOrdinal ) const
        {
            switch ( aEdgeOrdinal )
            {
                case ( 0 ):
                {
                    return { { 0, 1 } };
                    break;
                }
                case ( 1 ):
                {
                    return { { 1, 2 } };
                    break;
                }
                case ( 2 ):
                {
                    return { { 2, 3 } };
                    break;
                }
                case ( 3 ):
                {
                    return { { 3, 0 } };
                    break;
                }
                case ( 4 ):
                {
                    return { { 4, 5 } };
                    break;
                }
                case ( 5 ):
                {
                    return { { 5, 6 } };
                    break;
                }
                case ( 6 ):
                {
                    return { { 6, 7 } };
                    break;
                }
                case ( 7 ):
                {
                    return { { 7, 4 } };
                    break;
                }
                case ( 8 ):
                {
                    return { { 0, 4 } };
                    break;
                }
                case ( 9 ):
                {
                    return { { 1, 5 } };
                    break;
                }
                case ( 10 ):
                {
                    return { { 2, 6 } };
                    break;
                }
                case ( 11 ):
                {
                    return { { 3, 7 } };
                    break;
                }
                default:
                    MORIS_ASSERT( 0, "Invalid edge ordinal specified" );
                    return moris::Matrix< moris::IndexMat >( 0, 0 );
                    break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_to_facet_map( moris::uint aSideOrdinal ) const
        {
            return this->get_node_to_face_map( aSideOrdinal );
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_geometric_node_to_facet_map() const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_geometric_node_to_facet_map( moris::uint aSideOrdinal ) const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map( aSideOrdinal );
        }

        // ----------------------------------------------------------------------------------

        moris::uint
        Cell_Info_Hex20::get_adjacent_side_ordinal( moris::uint aSideOrdinal ) const
        {
            switch ( aSideOrdinal )
            {
                case 0:
                {
                    return 2;
                    break;
                }
                case 1:
                {
                    return 3;
                    break;
                }
                case 2:
                {
                    return 0;
                    break;
                }
                case 3:
                {
                    return 1;
                    break;
                }
                case 4:
                {
                    return 5;
                    break;
                }
                case 5:
                {
                    return 4;
                    break;
                }
                default:
                {
                    MORIS_ERROR( 0, "Invalid side ordinal specified" );
                    return MORIS_UINT_MAX;
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Vector< moris_index >
        Cell_Info_Hex20::get_vertex_path_to_entity_rank_and_ordinal(
                moris_index aVertexOrdinal,
                moris_index aOtherEntityOrdinal,
                moris_index aOtherEntityRank ) const
        {
            switch ( aOtherEntityRank )
            {
                // node to node paths
                case 0:
                {
                    Matrix< IndexMat > tVertexToVertexRanks = {
                        { -1, +1, +2, +1, +1, +2, +3, +2, +1, +2, +2, +1, +1, +2, +3, +2, +2, +3, +3, +2 },
                        { +1, -1, +1, +2, +2, +1, +2, +3, +1, +1, +2, +2, +2, +1, +2, +3, +2, +2, +3, +3 },
                        { +2, +1, -1, +1, +3, +2, +1, +2, +2, +1, +1, +2, +3, +2, +1, +2, +3, +2, +2, +3 },
                        { +1, +2, +1, -1, +2, +3, +2, +1, +2, +2, +1, +1, +2, +3, +2, +1, +3, +3, +2, +2 },
                        { +1, +2, +3, +2, -1, +1, +2, +1, +2, +3, +3, +2, +1, +2, +3, +2, +1, +2, +2, +1 },
                        { +2, +1, +2, +3, +1, -1, +1, +2, +2, +2, +3, +3, +2, +1, +2, +3, +1, +1, +2, +2 },
                        { +3, +2, +1, +2, +2, +1, -1, +1, +3, +2, +2, +3, +3, +2, +1, +2, +2, +1, +1, +2 },
                        { +2, +3, +2, +1, +1, +2, +1, -1, +3, +3, +2, +2, +2, +3, +2, +1, +2, +2, +1, +1 },
                        { +1, +1, +2, +2, +2, +2, +3, +3, -1, +2, +2, +2, +2, +2, +3, +3, +2, +3, +3, +3 },
                        { +2, +1, +1, +2, +3, +2, +2, +3, +2, -1, +2, +2, +3, +2, +2, +3, +3, +2, +3, +3 },
                        { +2, +2, +1, +1, +3, +3, +2, +2, +2, +2, -1, +2, +3, +3, +2, +2, +3, +3, +2, +3 },
                        { +1, +2, +2, +1, +2, +3, +3, +2, +2, +2, +2, -1, +2, +3, +3, +2, +3, +3, +3, +2 },
                        { +1, +2, +3, +2, +1, +2, +3, +2, +2, +3, +3, +2, -1, +2, +3, +2, +2, +3, +3, +2 },
                        { +2, +1, +2, +3, +2, +1, +2, +3, +2, +2, +3, +3, +2, -1, +2, +3, +2, +2, +3, +3 },
                        { +3, +2, +1, +2, +3, +2, +1, +2, +3, +2, +2, +3, +3, +2, -1, +2, +3, +2, +2, +3 },
                        { +2, +3, +2, +1, +2, +3, +2, +1, +3, +3, +2, +2, +2, +3, +2, -1, +3, +3, +2, +2 },
                        { +2, +2, +3, +3, +1, +1, +2, +2, +2, +3, +3, +3, +2, +2, +3, +3, -1, +2, +2, +2 },
                        { +3, +2, +2, +3, +2, +1, +1, +2, +3, +2, +3, +3, +3, +2, +2, +3, +2, -1, +2, +2 },
                        { +3, +3, +2, +2, +2, +2, +1, +1, +3, +3, +2, +3, +3, +3, +2, +2, +2, +2, -1, +2 },
                        { +2, +3, +3, +2, +1, +2, +2, +1, +3, +3, +3, +2, +2, +3, +3, +2, +2, +2, +2, -1 }
                    };

                    Matrix< IndexMat > tVertexToVertexIndices = {
                        { -1, +0, +4, +3, +8, +0, +0, +3, +0, +4, +4, +3, +8, +0, +0, +3, +0, +0, +0, +3 },
                        { +0, -1, +1, +4, +0, +9, +1, +0, +0, +1, +4, +4, +0, +9, +1, +0, +0, +1, +0, +0 },
                        { +4, +1, -1, +2, +0, +1, +10, +2, +4, +1, +2, +4, +0, +1, +10, +2, +0, +1, +2, +0 },
                        { +3, +4, +2, -1, +3, +0, +2, +11, +4, +4, +2, +3, +3, +0, +2, +11, +0, +0, +2, +3 },
                        { +8, +0, +0, +3, -1, +4, +5, +7, +0, +0, +0, +3, +8, +0, +0, +3, +4, +5, +5, +7 },
                        { +0, +9, +1, +0, +4, -1, +5, +5, +0, +1, +0, +0, +0, +9, +1, +0, +4, +5, +5, +5 },
                        { +0, +1, +10, +2, +5, +5, -1, +6, +0, +1, +2, +0, +0, +1, +10, +2, +5, +5, +6, +5 },
                        { +3, +0, +2, +11, +7, +5, +6, -1, +0, +0, +2, +3, +3, +0, +2, +11, +5, +5, +6, +7 },
                        { +0, +0, +4, +4, +0, +0, +0, +0, -1, +4, +4, +4, +0, +0, +0, +0, +0, +0, +0, +0 },
                        { +4, +1, +1, +4, +0, +1, +1, +0, +4, -1, +4, +4, +0, +1, +1, +0, +0, +1, +0, +0 },
                        { +4, +4, +2, +2, +0, +0, +2, +2, +4, +4, -1, +4, +0, +0, +2, +2, +0, +0, +2, +0 },
                        { +3, +4, +4, +3, +3, +0, +0, +3, +4, +4, +4, -1, +3, +0, +0, +3, +0, +0, +0, +3 },
                        { +8, +0, +0, +3, +8, +0, +0, +3, +0, +0, +0, +3, -1, +0, +0, +3, +0, +0, +0, +3 },
                        { +0, +9, +1, +0, +0, +9, +1, +0, +0, +1, +0, +0, +0, -1, +1, +0, +0, +1, +0, +0 },
                        { +0, +1, +10, +2, +0, +1, +10, +2, +0, +1, +2, +0, +0, +1, -1, +2, +0, +1, +2, +0 },
                        { +3, +0, +2, +11, +3, +0, +2, +11, +0, +0, +2, +3, +3, +0, +2, -1, +0, +0, +2, +3 },
                        { +0, +0, +0, +0, +4, +4, +5, +5, +0, +0, +0, +0, +0, +0, +0, +0, -1, +5, +5, +5 },
                        { +0, +1, +1, +0, +5, +5, +5, +5, +0, +1, +0, +0, +0, +1, +1, +0, +5, -1, +5, +5 },
                        { +0, +0, +2, +2, +5, +5, +6, +6, +0, +0, +2, +0, +0, +0, +2, +2, +5, +5, -1, +5 },
                        { +3, +0, +0, +3, +7, +5, +5, +7, +0, +0, +0, +3, +3, +0, +0, +3, +5, +5, +5, -1 }
                    };

                    moris_index tPathRank  = tVertexToVertexRanks( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );
                    moris_index tPathIndex = tVertexToVertexIndices( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );

                    MORIS_ASSERT( tPathRank != -1 && tPathIndex != -1,
                            "Cell_Info_Hex20::get_vertex_path_to_entity_rank_and_ordinal() - Vertex doesn't have path to itself." );

                    return { tPathIndex, tPathRank };
                    break;
                }

                // node to edge paths
                case 1:
                {
                    Matrix< IndexMat > tVertexToEdgeRanks = {
                        { +1, +2, +2, +1, +2, +3, +3, +2, +1, +2, +3, +2 },
                        { +1, +1, +2, +2, +2, +2, +3, +3, +2, +1, +2, +3 },
                        { +2, +1, +1, +2, +3, +2, +2, +3, +3, +2, +1, +2 },
                        { +2, +2, +1, +1, +3, +3, +2, +2, +2, +3, +2, +1 },
                        { +2, +3, +3, +2, +1, +2, +2, +1, +1, +2, +3, +2 },
                        { +2, +2, +3, +3, +1, +1, +2, +2, +2, +1, +2, +3 },
                        { +3, +2, +2, +3, +2, +1, +1, +2, +3, +2, +1, +2 },
                        { +3, +3, +2, +2, +2, +2, +1, +1, +2, +3, +2, +1 },
                        { +1, +2, +2, +2, +2, +3, +3, +3, +2, +2, +3, +3 },
                        { +2, +1, +2, +2, +3, +2, +3, +3, +3, +2, +2, +3 },
                        { +2, +2, +1, +2, +3, +3, +2, +3, +3, +3, +2, +2 },
                        { +2, +2, +2, +1, +3, +3, +3, +2, +2, +3, +3, +2 },
                        { +2, +3, +3, +2, +2, +3, +3, +2, +1, +2, +3, +2 },
                        { +2, +2, +3, +3, +2, +2, +3, +3, +2, +1, +2, +3 },
                        { +3, +2, +2, +3, +3, +2, +2, +3, +3, +2, +1, +2 },
                        { +3, +3, +2, +2, +3, +3, +2, +2, +2, +3, +2, +1 },
                        { +2, +3, +3, +3, +1, +2, +2, +2, +2, +2, +3, +3 },
                        { +3, +2, +3, +3, +2, +1, +2, +2, +3, +2, +2, +3 },
                        { +3, +3, +2, +3, +2, +2, +1, +2, +3, +3, +2, +2 },
                        { +3, +3, +3, +2, +2, +2, +2, +1, +2, +3, +3, +2 }
                    };

                    Matrix< IndexMat > tVertexToEdgeIndices = {
                        { +0, +4, +4, +3, +0, +0, +0, +3, +8, +0, +0, +3 },
                        { +0, +1, +4, +4, +0, +1, +0, +0, +0, +9, +1, +0 },
                        { +4, +1, +2, +4, +0, +1, +2, +0, +0, +1, +10, +2 },
                        { +4, +4, +2, +3, +0, +0, +2, +3, +3, +0, +2, +11 },
                        { +0, +0, +0, +3, +4, +5, +5, +7, +8, +0, +0, +3 },
                        { +0, +1, +0, +0, +4, +5, +5, +5, +0, +9, +1, +0 },
                        { +0, +1, +2, +0, +5, +5, +6, +5, +0, +1, +10, +2 },
                        { +0, +0, +2, +3, +5, +5, +6, +7, +3, +0, +2, +11 },
                        { +0, +4, +4, +4, +0, +0, +0, +0, +0, +0, +0, +0 },
                        { +4, +1, +4, +4, +0, +1, +0, +0, +0, +1, +1, +0 },
                        { +4, +4, +2, +4, +0, +0, +2, +0, +0, +0, +2, +2 },
                        { +4, +4, +4, +3, +0, +0, +0, +3, +3, +0, +0, +3 },
                        { +0, +0, +0, +3, +0, +0, +0, +3, +8, +0, +0, +3 },
                        { +0, +1, +0, +0, +0, +1, +0, +0, +0, +9, +1, +0 },
                        { +0, +1, +2, +0, +0, +1, +2, +0, +0, +1, +10, +2 },
                        { +0, +0, +2, +3, +0, +0, +2, +3, +3, +0, +2, +11 },
                        { +0, +0, +0, +0, +4, +5, +5, +5, +0, +0, +0, +0 },
                        { +0, +1, +0, +0, +5, +5, +5, +5, +0, +1, +1, +0 },
                        { +0, +0, +2, +0, +5, +5, +6, +5, +0, +0, +2, +2 },
                        { +0, +0, +0, +3, +5, +5, +5, +7, +3, +0, +0, +3 }
                    };

                    moris_index tPathRank  = tVertexToEdgeRanks( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );
                    moris_index tPathIndex = tVertexToEdgeIndices( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );

                    MORIS_ASSERT( tPathRank != -1 && tPathIndex != -1,
                            "Cell_Info_Hex20::get_vertex_path_to_entity_rank_and_ordinal() - Vertex doesn't have path to itself." );

                    return { tPathIndex, tPathRank };
                    break;
                }

                // node to face paths
                case 2:
                {
                    Matrix< IndexMat > tVertexToFaceRanks = {
                        { +2, +3, +3, +2, +2, +3 },
                        { +2, +2, +3, +3, +2, +3 },
                        { +3, +2, +2, +3, +2, +3 },
                        { +3, +3, +2, +2, +2, +3 },
                        { +2, +3, +3, +2, +3, +2 },
                        { +2, +2, +3, +3, +3, +2 },
                        { +3, +2, +2, +3, +3, +2 },
                        { +3, +3, +2, +2, +3, +2 },
                        { +2, +3, +3, +3, +2, +3 },
                        { +3, +2, +3, +3, +2, +3 },
                        { +3, +3, +2, +3, +2, +3 },
                        { +3, +3, +3, +2, +2, +3 },
                        { +2, +3, +3, +2, +3, +3 },
                        { +2, +2, +3, +3, +3, +3 },
                        { +3, +2, +2, +3, +3, +3 },
                        { +3, +3, +2, +2, +3, +3 },
                        { +2, +3, +3, +3, +3, +2 },
                        { +3, +2, +3, +3, +3, +2 },
                        { +3, +3, +2, +3, +3, +2 },
                        { +3, +3, +3, +2, +3, +2 }
                    };

                    Matrix< IndexMat > tVertexToFaceIndices = {
                        { +0, +0, +0, +3, +4, +0 },
                        { +0, +1, +0, +0, +4, +0 },
                        { +0, +1, +2, +0, +4, +0 },
                        { +0, +0, +2, +3, +4, +0 },
                        { +0, +0, +0, +3, +0, +5 },
                        { +0, +1, +0, +0, +0, +5 },
                        { +0, +1, +2, +0, +0, +5 },
                        { +0, +0, +2, +3, +0, +5 },
                        { +0, +0, +0, +0, +4, +0 },
                        { +0, +1, +0, +0, +4, +0 },
                        { +0, +0, +2, +0, +4, +0 },
                        { +0, +0, +0, +3, +4, +0 },
                        { +0, +0, +0, +3, +0, +0 },
                        { +0, +1, +0, +0, +0, +0 },
                        { +0, +1, +2, +0, +0, +0 },
                        { +0, +0, +2, +3, +0, +0 },
                        { +0, +0, +0, +0, +0, +5 },
                        { +0, +1, +0, +0, +0, +5 },
                        { +0, +0, +2, +0, +0, +5 },
                        { +0, +0, +0, +3, +0, +5 }
                    };

                    moris_index tPathRank  = tVertexToFaceRanks( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );
                    moris_index tPathIndex = tVertexToFaceIndices( (uint)aVertexOrdinal, (uint)aOtherEntityOrdinal );

                    MORIS_ASSERT( tPathRank != -1 && tPathIndex != -1,
                            "Cell_Info_Hex20::get_vertex_path_to_entity_rank_and_ordinal() - Vertex doesn't have path to itself." );

                    return { tPathIndex, tPathRank };
                    break;
                }

                default:
                {
                    MORIS_ERROR( 0, "Invalid other entity rank for hex20" );
                    return moris::Vector< moris_index >( 0 );
                }
            }    // end: switch aOtherEntityRank
        }

        // ----------------------------------------------------------------------------------

        moris::Vector< moris_index >
        Cell_Info_Hex20::get_edge_path_to_entity_rank_and_ordinal(
                moris_index aEdgeOrdinal,
                moris_index aOtherEntityOrdinal,
                moris_index aOtherEntityRank ) const
        {
            switch ( aEdgeOrdinal )
            {
                case 0:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 1 };
                                case 1:
                                    return { 4, 2 };
                                case 2:
                                    return { 4, 2 };
                                case 3:
                                    return { 4, 2 };
                                case 4:
                                    return { 0, 2 };
                                case 5:
                                    return { 0, 3 };
                                case 6:
                                    return { 0, 3 };
                                case 7:
                                    return { 0, 3 };
                                case 8:
                                    return { 0, 2 };
                                case 9:
                                    return { 0, 2 };
                                case 10:
                                    return { 0, 3 };
                                case 11:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 4, 2 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 1:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 4, 2 };
                                case 1:
                                    return { 1, 1 };
                                case 2:
                                    return { 4, 2 };
                                case 3:
                                    return { 4, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 1, 2 };
                                case 6:
                                    return { 0, 3 };
                                case 7:
                                    return { 0, 3 };
                                case 8:
                                    return { 0, 3 };
                                case 9:
                                    return { 1, 2 };
                                case 10:
                                    return { 1, 2 };
                                case 11:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // edge to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 4, 2 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 2:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 4, 2 };
                                case 1:
                                    return { 4, 2 };
                                case 2:
                                    return { 2, 1 };
                                case 3:
                                    return { 4, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                case 6:
                                    return { 2, 2 };
                                case 7:
                                    return { 0, 3 };
                                case 8:
                                    return { 0, 3 };
                                case 9:
                                    return { 0, 3 };
                                case 10:
                                    return { 2, 2 };
                                case 11:
                                    return { 2, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 4, 2 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 3:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 4, 2 };
                                case 1:
                                    return { 4, 2 };
                                case 2:
                                    return { 4, 2 };
                                case 3:
                                    return { 3, 1 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                case 6:
                                    return { 0, 3 };
                                case 7:
                                    return { 3, 2 };
                                case 8:
                                    return { 3, 2 };
                                case 9:
                                    return { 0, 3 };
                                case 10:
                                    return { 0, 3 };
                                case 11:
                                    return { 3, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 4, 2 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 4:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 4, 1 };
                                case 5:
                                    return { 5, 2 };
                                case 6:
                                    return { 5, 2 };
                                case 7:
                                    return { 5, 2 };
                                case 8:
                                    return { 0, 2 };
                                case 9:
                                    return { 0, 2 };
                                case 10:
                                    return { 0, 3 };
                                case 11:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 5, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 5:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 5, 2 };
                                case 5:
                                    return { 5, 1 };
                                case 6:
                                    return { 5, 2 };
                                case 7:
                                    return { 5, 2 };
                                case 8:
                                    return { 0, 3 };
                                case 9:
                                    return { 1, 2 };
                                case 10:
                                    return { 1, 2 };
                                case 11:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 5, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 6:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 5, 2 };
                                case 5:
                                    return { 5, 2 };
                                case 6:
                                    return { 6, 1 };
                                case 7:
                                    return { 5, 2 };
                                case 8:
                                    return { 0, 3 };
                                case 9:
                                    return { 0, 3 };
                                case 10:
                                    return { 2, 2 };
                                case 11:
                                    return { 2, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 5, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 7:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 5, 2 };
                                case 5:
                                    return { 5, 2 };
                                case 6:
                                    return { 5, 2 };
                                case 7:
                                    return { 7, 1 };
                                case 8:
                                    return { 3, 2 };
                                case 9:
                                    return { 0, 3 };
                                case 10:
                                    return { 0, 3 };
                                case 11:
                                    return { 3, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 5, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 8:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 0, 2 };
                                case 5:
                                    return { 0, 3 };
                                case 6:
                                    return { 0, 3 };
                                case 7:
                                    return { 3, 2 };
                                case 8:
                                    return { 8, 1 };
                                case 9:
                                    return { 0, 2 };
                                case 10:
                                    return { 0, 3 };
                                case 11:
                                    return { 3, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 9:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 2 };
                                case 5:
                                    return { 1, 2 };
                                case 6:
                                    return { 0, 3 };
                                case 7:
                                    return { 0, 3 };
                                case 8:
                                    return { 0, 2 };
                                case 9:
                                    return { 9, 1 };
                                case 10:
                                    return { 1, 2 };
                                case 11:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 2 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 0, 3 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 10:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 1, 2 };
                                case 6:
                                    return { 2, 2 };
                                case 7:
                                    return { 0, 3 };
                                case 8:
                                    return { 0, 3 };
                                case 9:
                                    return { 1, 2 };
                                case 10:
                                    return { 10, 1 };
                                case 11:
                                    return { 2, 2 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 1, 2 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 0, 3 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                case 11:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                case 6:
                                    return { 2, 2 };
                                case 7:
                                    return { 3, 2 };
                                case 8:
                                    return { 3, 2 };
                                case 9:
                                    return { 0, 3 };
                                case 10:
                                    return { 2, 2 };
                                case 11:
                                    return { 11, 1 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }

                        // node to face paths
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return { 0, 3 };
                                case 1:
                                    return { 0, 3 };
                                case 2:
                                    return { 2, 2 };
                                case 3:
                                    return { 3, 2 };
                                case 4:
                                    return { 0, 3 };
                                case 5:
                                    return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Vector< moris_index >( 0 );
                                }
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Vector< moris_index >( 0 );
                        }
                    }
                }
                default:
                {
                    MORIS_ERROR( 0, "Invalid vertex ordinal for hex8" );
                    return moris::Vector< moris_index >( 0 );
                }
            }
        }

        // ----------------------------------------------------------------------------------

        bool
        Cell_Info_Hex20::is_entity_connected_to_facet(
                moris_index aFacetOrdinal,
                moris_index aOtherEntityOrdinal,
                moris_index aOtherEntityRank ) const
        {
            switch ( aFacetOrdinal )
            {
                case 0:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                case 1:
                                    return true;
                                case 2:
                                    return false;
                                case 3:
                                    return false;
                                case 4:
                                    return true;
                                case 5:
                                    return true;
                                case 6:
                                    return false;
                                case 7:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                case 1:
                                    return false;
                                case 2:
                                    return false;
                                case 3:
                                    return false;
                                case 4:
                                    return true;
                                case 5:
                                    return false;
                                case 6:
                                    return false;
                                case 7:
                                    return false;
                                case 8:
                                    return true;
                                case 9:
                                    return true;
                                case 10:
                                    return false;
                                case 11:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex20" );
                            return false;
                        }
                    }
                }
                case 1:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return true;
                                case 2:
                                    return true;
                                case 3:
                                    return false;
                                case 4:
                                    return false;
                                case 5:
                                    return true;
                                case 6:
                                    return true;
                                case 7:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return true;
                                case 2:
                                    return false;
                                case 3:
                                    return false;
                                case 4:
                                    return false;
                                case 5:
                                    return true;
                                case 6:
                                    return false;
                                case 7:
                                    return false;
                                case 8:
                                    return false;
                                case 9:
                                    return true;
                                case 10:
                                    return true;
                                case 11:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 1:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex20" );
                            return false;
                        }
                    }
                }
                case 2:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return false;
                                case 2:
                                    return true;
                                case 3:
                                    return true;
                                case 4:
                                    return false;
                                case 5:
                                    return false;
                                case 6:
                                    return true;
                                case 7:
                                    return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return false;
                                case 2:
                                    return true;
                                case 3:
                                    return false;
                                case 4:
                                    return false;
                                case 5:
                                    return false;
                                case 6:
                                    return true;
                                case 7:
                                    return false;
                                case 8:
                                    return false;
                                case 9:
                                    return false;
                                case 10:
                                    return true;
                                case 11:
                                    return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 2:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return false;
                        }
                    }
                }
                case 3:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                case 1:
                                    return false;
                                case 2:
                                    return false;
                                case 3:
                                    return true;
                                case 4:
                                    return true;
                                case 5:
                                    return false;
                                case 6:
                                    return false;
                                case 7:
                                    return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return false;
                                case 2:
                                    return false;
                                case 3:
                                    return true;
                                case 4:
                                    return false;
                                case 5:
                                    return false;
                                case 6:
                                    return false;
                                case 7:
                                    return true;
                                case 8:
                                    return true;
                                case 9:
                                    return false;
                                case 10:
                                    return false;
                                case 11:
                                    return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex20" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 3:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex20" );
                            return false;
                        }
                    }
                }
                case 4:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                case 1:
                                    return true;
                                case 2:
                                    return true;
                                case 3:
                                    return true;
                                case 4:
                                    return false;
                                case 5:
                                    return false;
                                case 6:
                                    return false;
                                case 7:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return true;
                                case 1:
                                    return true;
                                case 2:
                                    return true;
                                case 3:
                                    return true;
                                case 4:
                                    return false;
                                case 5:
                                    return false;
                                case 6:
                                    return false;
                                case 7:
                                    return false;
                                case 8:
                                    return false;
                                case 9:
                                    return false;
                                case 10:
                                    return false;
                                case 11:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 4:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return false;
                        }
                    }
                }
                case 5:
                {
                    switch ( aOtherEntityRank )
                    {
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return false;
                                case 2:
                                    return false;
                                case 3:
                                    return false;
                                case 4:
                                    return true;
                                case 5:
                                    return true;
                                case 6:
                                    return true;
                                case 7:
                                    return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return false;
                                }
                            }
                        }
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0:
                                    return false;
                                case 1:
                                    return false;
                                case 2:
                                    return false;
                                case 3:
                                    return false;
                                case 4:
                                    return true;
                                case 5:
                                    return true;
                                case 6:
                                    return true;
                                case 7:
                                    return true;
                                case 8:
                                    return false;
                                case 9:
                                    return false;
                                case 10:
                                    return false;
                                case 11:
                                    return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return false;
                                }
                            }
                        }
                        case 2:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 5:
                                    return true;
                                default:
                                    return false;
                            }
                        }
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return false;
                        }
                    }
                }
                default:
                {
                    MORIS_ERROR( 0, "Invalid facet ordinal for hex8" );
                    return false;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        Matrix< DDRMat >
        Cell_Info_Hex20::get_vertex_loc_coord( moris_index const & aVertexOrdinal ) const
        {
            switch ( aVertexOrdinal )
            {
                case 0:
                {
                    return { { -1.0, -1.0, -1.0 } };
                    break;
                }
                case 1:
                {
                    return { { +1.0, -1.0, -1.0 } };
                    break;
                }
                case 2:
                {
                    return { { +1.0, +1.0, -1.0 } };
                    break;
                }
                case 3:
                {
                    return { { -1.0, +1.0, -1.0 } };
                    break;
                }
                case 4:
                {
                    return { { -1.0, -1.0, +1.0 } };
                    break;
                }
                case 5:
                {
                    return { { +1.0, -1.0, +1.0 } };
                    break;
                }
                case 6:
                {
                    return { { +1.0, +1.0, +1.0 } };
                    break;
                }
                case 7:
                {
                    return { { -1.0, +1.0, +1.0 } };
                    break;
                }
                case 8:
                {
                    return { { +0.0, -1.0, -1.0 } };
                    break;
                }
                case 9:
                {
                    return { { +1.0, 0.0, -1.0 } };
                    break;
                }
                case 10:
                {
                    return { { 0.0, +1.0, -1.0 } };
                    break;
                }
                case 11:
                {
                    return { { -1.0, 0.0, -1.0 } };
                    break;
                }
                case 12:
                {
                    return { { -1.0, -1.0, 0.0 } };
                    break;
                }
                case 13:
                {
                    return { { +1.0, -1.0, 0.0 } };
                    break;
                }
                case 14:
                {
                    return { { +1.0, +1.0, 0.0 } };
                    break;
                }
                case 15:
                {
                    return { { -1.0, +1.0, 0.0 } };
                    break;
                }
                case 16:
                {
                    return { { 0.0, -1.0, +1.0 } };
                    break;
                }
                case 17:
                {
                    return { { +1.0, 0.0, +1.0 } };
                    break;
                }
                case 18:
                {
                    return { { 0.0, +1.0, +1.0 } };
                    break;
                }
                case 19:
                {
                    return { { -1.0, 0.0, +1.0 } };
                    break;
                }
                default:
                {
                    MORIS_ERROR( 0, "Invalid vertex ordinal specified" );
                    return moris::Matrix< moris::DDRMat >( 0, 0 );
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------
        inline moris::Matrix< moris::IndexMat >
        Cell_Info_Hex20::get_node_map_outward_normal( moris::uint aSideOrdinal ) const
        {
            switch ( aSideOrdinal )
            {
                case ( 0 ):
                {
                    return { { 0, 1 }, { 1, 5 } };
                    break;
                }
                case ( 1 ):
                {
                    return { { 1, 2 }, { 2, 6 } };
                    break;
                }
                case ( 2 ):
                {
                    return { { 2, 3 }, { 3, 7 } };
                    break;
                }
                case ( 3 ):
                {
                    return { { 7, 3 }, { 3, 0 } };
                    break;
                }
                case ( 4 ):
                {
                    return { { 0, 3 }, { 3, 2 } };
                    break;
                }
                case ( 5 ):
                {
                    return { { 4, 5 }, { 5, 6 } };
                    break;
                }
                default:
                    MORIS_ERROR( 0, "Invalid side ordinal specified" );
                    return moris::Matrix< moris::IndexMat >( 0, 0 );
                    break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex20::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            moris::Vector< Vertex* > tVertices = aCell->get_vertex_pointers();

            // FIXME: only works for rectangular cells
            Matrix< DDRMat > tNode0Coords = tVertices( 0 )->get_coords();
            Matrix< DDRMat > tNode6Coords = tVertices( 6 )->get_coords();

            real tLx = std::abs( tNode0Coords( 0 ) - tNode6Coords( 0 ) );
            real tLy = std::abs( tNode0Coords( 1 ) - tNode6Coords( 1 ) );
            real tLz = std::abs( tNode0Coords( 2 ) - tNode6Coords( 2 ) );

            return tLx * tLy * tLz;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex20::compute_cell_side_size(
                moris::mtk::Cell const * aCell,
                moris_index const &      aSideOrd ) const
        {
            moris::Vector< Vertex const * > tVertices = aCell->get_vertices_on_side_ordinal( aSideOrd );

            // FIXME: only works for rectangular cells
            Matrix< DDRMat > tNodeCoords0 = tVertices( 0 )->get_coords();
            Matrix< DDRMat > tNodeCoords1 = tVertices( 1 )->get_coords();
            Matrix< DDRMat > tNodeCoords2 = tVertices( 3 )->get_coords();

            return norm( cross( tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0 ) );
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Hex20::eval_N(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX20 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real xi   = aXi( 0 );
            real eta  = aXi( 1 );
            real zeta = aXi( 2 );

            // often used constants
            real xi2   = std::pow( xi, 2 );
            real eta2  = std::pow( eta, 2 );
            real zeta2 = std::pow( zeta, 2 );

            // populate output matrix
            aNXi.set_size( 1, 20 );
            aNXi( 0 )  = ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * ( eta + xi + zeta + 2.0 ) * 0.125;
            aNXi( 1 )  = -( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta - xi + zeta + 2.0 ) * 0.125;
            aNXi( 2 )  = -( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * ( eta + xi - zeta - 2.0 ) * 0.125;
            aNXi( 3 )  = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * ( -eta + xi + zeta + 2.0 ) * 0.125;
            aNXi( 4 )  = -( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta + xi - zeta + 2.0 ) * 0.125;
            aNXi( 5 )  = ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta - xi - zeta + 2.0 ) * 0.125;
            aNXi( 6 )  = ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * ( eta + xi + zeta - 2.0 ) * 0.125;
            aNXi( 7 )  = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * ( eta - xi + zeta - 2.0 ) * 0.125;
            aNXi( 8 )  = -( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 9 )  = ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 10 ) = ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 11 ) = -( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.25;
            aNXi( 12 ) = -( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 ) * 0.25;
            aNXi( 13 ) = ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 ) * 0.25;
            aNXi( 14 ) = -( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 ) * 0.25;
            aNXi( 15 ) = ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 ) * 0.25;
            aNXi( 16 ) = ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 17 ) = -( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 18 ) = -( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 ) * 0.25;
            aNXi( 19 ) = ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.25;
        }

    }    // namespace mtk
}    // namespace moris
