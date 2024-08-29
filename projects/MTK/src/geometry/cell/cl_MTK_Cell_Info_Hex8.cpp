/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Hex8.cpp
 *
 */

#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "cl_Vector.hpp"

namespace moris::mtk
{
    // ----------------------------------------------------------------------------------

    enum Geometry_Type
    Cell_Info_Hex8::get_cell_geometry() const
    {
        return Geometry_Type::HEX;
    }

    // ----------------------------------------------------------------------------------

    enum CellTopology
    Cell_Info_Hex8::get_cell_topology() const
    {
        return CellTopology::HEX8;
    }

    // ----------------------------------------------------------------------------------

    enum Interpolation_Order
    Cell_Info_Hex8::get_cell_interpolation_order() const
    {
        return Interpolation_Order::LINEAR;
    }

    // ----------------------------------------------------------------------------------

    enum Integration_Order
    Cell_Info_Hex8::get_cell_integration_order() const
    {
        return Integration_Order::HEX_2x2x2;
    }

    //-----------------------------------------------------------------------------

    enum CellShape
    Cell_Info_Hex8::compute_cell_shape( moris::mtk::Cell const *aCell ) const
    {
        // getting vertices and storing them in a local matrix, since each node will be used a few times
        Vector< Vertex * > tVertices = aCell->get_vertex_pointers();

        // init cell shape
        CellShape tCellShape = CellShape::RECTANGULAR;

        // error threshold
        real tEpsilon = 1.0E-8;

        // init cell of face normals
        Vector< moris::Matrix< DDRMat > > tFaceNormals( 6 );

        // looping through each face
        for ( uint iFace = 0; iFace < 6; iFace++ )
        {
            // getting nodes on the face
            moris::Matrix< moris::IndexMat > tFaceNodes = this->get_node_to_face_map( iFace );

            // get all of the face vertices
            moris::Matrix< DDRMat > tVertex0 = tVertices( tFaceNodes( 0 ) )->get_coords();
            moris::Matrix< DDRMat > tVertex1 = tVertices( tFaceNodes( 1 ) )->get_coords();
            moris::Matrix< DDRMat > tVertex2 = tVertices( tFaceNodes( 2 ) )->get_coords();
            moris::Matrix< DDRMat > tVertex3 = tVertices( tFaceNodes( 3 ) )->get_coords();

            // get edges to define check plane
            auto tEdge0           = tVertex1 - tVertex0;
            auto tEdge1           = tVertex2 - tVertex1;
            tFaceNormals( iFace ) = cross( tEdge0, tEdge1 );

            // these are the other plane normals produced by the other 6 nodes.
            // keeping node 0 to ensure the plane is at the same offset.
            auto tEdge2     = tVertex1 - tVertex0;
            auto tEdge3     = tVertex3 - tVertex1;
            auto tFaceVec23 = cross( tEdge2, tEdge3 );

            // All three of the plane normals must be parallel in order to be considered straight shape
            if ( norm( cross( tFaceNormals( iFace ), tFaceVec23 ) ) > tEpsilon )
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
            if ( norm( tCross02 ) > tEpsilon || norm( tCross13 ) > tEpsilon || norm( tCross45 ) > tEpsilon )
            {
                tCellShape = CellShape::STRAIGHT;
            }
            else
            {
                // get the some edges
                moris::Matrix< DDRMat > tVertex1 = tVertices( 1 )->get_coords();
                moris::Matrix< DDRMat > tEdge0   = tVertex1 - tVertices( 0 )->get_coords();
                moris::Matrix< DDRMat > tEdge1   = tVertices( 2 )->get_coords() - tVertex1;

                // checking that this cell fits requirements for being rectangular and aligned (CellShape::RECTANGULAR)
                if ( std::abs( dot( tFaceNormals( 0 ), tFaceNormals( 1 ) ) ) > tEpsilon || std::abs( dot( tFaceNormals( 0 ), tFaceNormals( 4 ) ) ) > tEpsilon || std::abs( dot( tFaceNormals( 1 ), tFaceNormals( 4 ) ) ) > tEpsilon || std::abs( tEdge0( 1 ) ) > tEpsilon || std::abs( tEdge0( 2 ) ) > tEpsilon || std::abs( tEdge1( 0 ) ) > tEpsilon || std::abs( tEdge1( 2 ) ) > tEpsilon || tEdge0( 0 ) < 0.0 || tEdge1( 1 ) < 0.0 )
                {
                    tCellShape = CellShape::PARALLEL;
                }
            }
        }

        return tCellShape;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cell_Info_Hex8::get_num_verts() const
    {
        return 8;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cell_Info_Hex8::get_num_facets() const
    {
        return 6;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cell_Info_Hex8::get_num_edges() const
    {
        return 12;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cell_Info_Hex8::get_num_verts_per_facet() const
    {
        return 4;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cell_Info_Hex8::get_loc_coord_dim() const
    {
        return 3;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_node_to_face_map() const
    {
        return {
            { 1, 5, 4, 0 },
            { 1, 2, 6, 5 },
            { 3, 7, 6, 2 },
            { 0, 4, 7, 3 },
            { 0, 3, 2, 1 },
            { 4, 5, 6, 7 }
        };
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_node_to_edge_map() const
    {
        return {
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
        };
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_node_to_face_map( moris::uint aSideOrdinal ) const
    {
        switch ( aSideOrdinal )
        {
            case ( 0 ):
            {
                return { { 1, 5, 4, 0 } };
                break;
            }
            case ( 1 ):
            {
                return { { 1, 2, 6, 5 } };
                break;
            }
            case ( 2 ):
            {
                return { { 3, 7, 6, 2 } };
                break;
            }
            case ( 3 ):
            {
                return { { 0, 4, 7, 3 } };
                break;
            }
            case ( 4 ):
            {
                return { { 0, 3, 2, 1 } };
                break;
            }
            case ( 5 ):
            {
                return { { 4, 5, 6, 7 } };
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
    Cell_Info_Hex8::get_node_to_edge_map( moris::uint aEdgeOrdinal ) const
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
    Cell_Info_Hex8::get_node_to_facet_map() const
    {
        return this->get_node_to_face_map();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_node_to_facet_map( moris::uint aSideOrdinal ) const
    {
        return this->get_node_to_face_map( aSideOrdinal );
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_geometric_node_to_facet_map() const
    {
        return this->get_node_to_face_map();
    }

    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cell_Info_Hex8::get_vertex_loc_coord( moris_index const &aVertexOrdinal ) const
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
            default:
            {
                MORIS_ERROR( 0, "Invalid vertex ordinal specified" );
                return moris::Matrix< moris::DDRMat >( 0, 0 );
                break;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_geometric_node_to_facet_map( moris::uint aSideOrdinal ) const
    {
        return this->get_node_to_face_map( aSideOrdinal );
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Cell_Info_Hex8::get_node_map_outward_normal( moris::uint aSideOrdinal ) const
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
    Cell_Info_Hex8::compute_cell_size_special( moris::mtk::Cell const *aCell ) const
    {
        Vector< Vertex * > tVertices = aCell->get_vertex_pointers();

        Matrix< DDRMat > tNode0Coords = tVertices( 0 )->get_coords();
        Matrix< DDRMat > tNode6Coords = tVertices( 6 )->get_coords();

        // FIXME: only works for rectangular cells
        real tLx = std::abs( tNode0Coords( 0 ) - tNode6Coords( 0 ) );
        real tLy = std::abs( tNode0Coords( 1 ) - tNode6Coords( 1 ) );
        real tLz = std::abs( tNode0Coords( 2 ) - tNode6Coords( 2 ) );

        return tLx * tLy * tLz;
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Cell_Info_Hex8::compute_cell_size_straight( moris::mtk::Cell const *aCell ) const
    {
        // FIXME: Not consistent with each vertex. Depending what corners are used to be bisected, this won't be consistent

        // getting vertices and storing them in a local matrix, since each node will be used a few times
        Vector< Vertex * > tVertices = aCell->get_vertex_pointers();
        Matrix< DDRMat >   tNodeCoords( 8, 3 );

        for ( uint i = 0; i < 8; ++i )
        {
            tNodeCoords( { i, i }, { 0, 2 } ) = tVertices( i )->get_coords()( { 0, 0 }, { 0, 2 } );
        }

        // permutation matrix of vertex IDs to stipulate individual tet4 calculations
        moris::Matrix< DDUMat > tPermutMap = { { 1, 0, 5, 2 },
            { 4, 0, 7, 5 },
            { 3, 0, 2, 7 },
            { 6, 2, 5, 7 },
            { 7, 0, 2, 5 } };

        // init volume and tet edge vectors
        moris::real      tVolume = 0.0;
        Matrix< DDRMat > tNodeCoords10( 1, 3 );
        Matrix< DDRMat > tNodeCoords20( 1, 3 );
        Matrix< DDRMat > tNodeCoords30( 1, 3 );

        // A hex can be broken into 5 separate tets
        for ( uint iTet = 0; iTet < 5; ++iTet )
        {
            // Assigning Vectors
            tNodeCoords10 = tNodeCoords( { tPermutMap( iTet, 1 ), tPermutMap( iTet, 1 ) }, { 0, 2 } ) - tNodeCoords( { tPermutMap( iTet, 0 ), tPermutMap( iTet, 0 ) }, { 0, 2 } );
            tNodeCoords20 = tNodeCoords( { tPermutMap( iTet, 2 ), tPermutMap( iTet, 2 ) }, { 0, 2 } ) - tNodeCoords( { tPermutMap( iTet, 0 ), tPermutMap( iTet, 0 ) }, { 0, 2 } );
            tNodeCoords30 = tNodeCoords( { tPermutMap( iTet, 3 ), tPermutMap( iTet, 3 ) }, { 0, 2 } ) - tNodeCoords( { tPermutMap( iTet, 0 ), tPermutMap( iTet, 0 ) }, { 0, 2 } );

            MORIS_ASSERT( 1.0 / 6.0 * dot( tNodeCoords10, cross( tNodeCoords20, tNodeCoords30 ) ) > 0,
                    "Cell_Info_Hex8::compute_cell_size_straight - Determined interior tet "
                    "volume is <=0, suggesting poorly defined nodal coordinates." );

            tVolume += 1.0 / 6.0 * dot( tNodeCoords10, cross( tNodeCoords20, tNodeCoords30 ) );
        }

        return tVolume;
    }

    // ----------------------------------------------------------------------------------
    moris::real
    Cell_Info_Hex8::compute_cell_side_size( moris::mtk::Cell const *aCell,
            moris_index const                                      &aSideOrd ) const
    {
        Vector< mtk::Vertex const * > tVertices = aCell->get_vertices_on_side_ordinal( aSideOrd );

        // FIXME: only works for rectangular cells
        Matrix< DDRMat > tNodeCoords0 = tVertices( 0 )->get_coords();
        Matrix< DDRMat > tNodeCoords1 = tVertices( 1 )->get_coords();
        Matrix< DDRMat > tNodeCoords2 = tVertices( 3 )->get_coords();

        return norm( cross( tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0 ) );
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Cell_Info_Hex8::get_adjacent_side_ordinal( moris::uint aSideOrdinal ) const
    {
        switch ( aSideOrdinal )
        {
            case ( 0 ):
            {
                return 2;
                break;
            }
            case ( 1 ):
            {
                return 3;
                break;
            }
            case ( 2 ):
            {
                return 0;
                break;
            }
            case ( 3 ):
            {
                return 1;
                break;
            }
            case ( 4 ):
            {
                return 5;
                break;
            }
            case ( 5 ):
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

    Vector< moris_index >
    Cell_Info_Hex8::get_vertex_path_to_entity_rank_and_ordinal(
            moris_index aVertexOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
    {
        switch ( aVertexOrdinal )
        {
            case 0:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 1:
                                return { 0, 1 };
                            case 2:
                                return { 4, 2 };
                            case 3:
                                return { 3, 1 };
                            case 4:
                                return { 8, 1 };
                            case 5:
                                return { 0, 2 };
                            case 6:
                                return { 0, 3 };
                            case 7:
                                return { 3, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 3, 1 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 4, 2 };
                            case 5:
                                return { 0, 3 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 1:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 0, 1 };
                            case 1:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 2:
                                return { 1, 1 };
                            case 3:
                                return { 4, 2 };
                            case 4:
                                return { 0, 2 };
                            case 5:
                                return { 9, 1 };
                            case 6:
                                return { 1, 2 };
                            case 7:
                                return { 0, 3 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
                    case 1:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 0, 1 };
                            case 1:
                                return { 1, 1 };
                            case 2:
                                return { 4, 2 };
                            case 3:
                                return { 4, 2 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 4, 2 };
                            case 5:
                                return { 0, 3 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 2:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 4, 2 };
                            case 1:
                                return { 1, 1 };
                            case 2:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 3:
                                return { 2, 1 };
                            case 4:
                                return { 0, 3 };
                            case 5:
                                return { 1, 2 };
                            case 6:
                                return { 10, 1 };
                            case 7:
                                return { 2, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
                    case 1:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 4, 2 };
                            case 1:
                                return { 1, 1 };
                            case 2:
                                return { 2, 1 };
                            case 3:
                                return { 4, 2 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 4, 2 };
                            case 5:
                                return { 0, 3 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 3:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 3, 1 };
                            case 1:
                                return { 4, 2 };
                            case 2:
                                return { 2, 1 };
                            case 3:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 4:
                                return { 3, 2 };
                            case 5:
                                return { 0, 3 };
                            case 6:
                                return { 2, 2 };
                            case 7:
                                return { 11, 1 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 3, 1 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 4, 2 };
                            case 5:
                                return { 0, 3 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 4:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 8, 1 };
                            case 1:
                                return { 0, 2 };
                            case 2:
                                return { 0, 3 };
                            case 3:
                                return { 3, 2 };
                            case 4:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 5:
                                return { 4, 1 };
                            case 6:
                                return { 5, 2 };
                            case 7:
                                return { 7, 1 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 4, 1 };
                            case 5:
                                return { 5, 2 };
                            case 6:
                                return { 5, 2 };
                            case 7:
                                return { 7, 1 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 5, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 5:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 0, 2 };
                            case 1:
                                return { 9, 1 };
                            case 2:
                                return { 1, 2 };
                            case 3:
                                return { 0, 3 };
                            case 4:
                                return { 4, 1 };
                            case 5:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 6:
                                return { 5, 1 };
                            case 7:
                                return { 5, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 4, 1 };
                            case 5:
                                return { 5, 1 };
                            case 6:
                                return { 5, 2 };
                            case 7:
                                return { 5, 2 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 5, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 6:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 0, 3 };
                            case 1:
                                return { 1, 2 };
                            case 2:
                                return { 10, 1 };
                            case 3:
                                return { 2, 2 };
                            case 4:
                                return { 5, 2 };
                            case 5:
                                return { 5, 1 };
                            case 6:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            case 7:
                                return { 6, 1 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 5, 2 };
                            case 5:
                                return { 5, 1 };
                            case 6:
                                return { 6, 1 };
                            case 7:
                                return { 5, 2 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 5, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            case 7:
            {
                switch ( aOtherEntityRank )
                {
                    // node to node paths
                    case 0:
                    {
                        switch ( aOtherEntityOrdinal )
                        {
                            case 0:
                                return { 3, 2 };
                            case 1:
                                return { 0, 3 };
                            case 2:
                                return { 2, 2 };
                            case 3:
                                return { 11, 1 };
                            case 4:
                                return { 7, 1 };
                            case 5:
                                return { 5, 2 };
                            case 6:
                                return { 6, 1 };
                            case 7:
                            {
                                MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                return 0;
                            }
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    // node to edge paths
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
                                return { 5, 2 };
                            case 5:
                                return { 5, 2 };
                            case 6:
                                return { 6, 1 };
                            case 7:
                                return { 7, 1 };
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
                                return Vector< moris_index >( 0 );
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
                                return { 5, 2 };
                            default:
                            {
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            default:
            {
                MORIS_ERROR( 0, "Invalid vertex ordinal for hex8" );
                return Vector< moris_index >( 0 );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    Vector< moris_index >
    Cell_Info_Hex8::get_edge_path_to_entity_rank_and_ordinal(
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
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
                                return Vector< moris_index >( 0 );
                            }
                        }
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                        return Vector< moris_index >( 0 );
                    }
                }
            }
            default:
            {
                MORIS_ERROR( 0, "Invalid vertex ordinal for hex8" );
                return Vector< moris_index >( 0 );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Cell_Info_Hex8::is_entity_connected_to_facet(
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
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
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
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
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
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
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
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
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
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
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
                                MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
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
                        MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
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

    void
    Cell_Info_Hex8::eval_N(
            const Matrix< DDRMat > &aXi,
            Matrix< DDRMat >       &aNXi ) const
    {
        // make sure that input is correct
        MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size." );

        // unpack xi and eta from input vector
        moris::real xi   = aXi( 0 );
        moris::real eta  = aXi( 1 );
        moris::real zeta = aXi( 2 );

        // populate output matrix
        aNXi.set_size( 1, 8 );
        aNXi( 0 ) = -( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aNXi( 1 ) = ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aNXi( 2 ) = -( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aNXi( 3 ) = ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aNXi( 4 ) = ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aNXi( 5 ) = -( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aNXi( 6 ) = ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aNXi( 7 ) = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
    }
}    // namespace moris::mtk
