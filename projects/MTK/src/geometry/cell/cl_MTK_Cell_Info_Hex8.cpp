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
        return this->compute_hex_cell_shape( aCell );
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
    Cell_Info_Hex8::compute_cell_size_deriv(
            moris::mtk::Cell const *aCell,
            uint                    aLocalVertexID,
            uint                    aDirection ) const
    {
        MORIS_ASSERT( aLocalVertexID < 8,
                "Cell_Info_Hex8::compute_cell_size_deriv - Vertex ordinal must be less than 8." );
        MORIS_ASSERT( aDirection < 3,
                "Cell_Info_Hex8::compute_cell_size_deriv - Direction must be 0, 1, or 2." );

        const Matrix< DDRMat > tCoordinates = aCell->get_vertex_coords();
        MORIS_ASSERT( tCoordinates.n_rows() == 8 && tCoordinates.n_cols() == 3,
                "Cell_Info_Hex8::compute_cell_size_deriv - Expected eight vertices in three dimensions." );

        // Signs of (xi, eta, zeta) in the standard HEX8 vertex ordering.
        constexpr real tSigns[ 8 ][ 3 ] = {
            { -1.0, -1.0, -1.0 }, { 1.0, -1.0, -1.0 }, { 1.0, 1.0, -1.0 }, { -1.0, 1.0, -1.0 }, { -1.0, -1.0, 1.0 }, { 1.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0 }, { -1.0, 1.0, 1.0 }
        };

        // Write x = a0 + a1*xi + a2*eta + a3*zeta + a4*xi*eta
        //             + a5*eta*zeta + a6*xi*zeta + a7*xi*eta*zeta.
        // Subtracting vertex 0 avoids cancellation from a large coordinate origin.
        // Only a1,...,a6 enter the integrated determinant; a0 and a7 cancel.
        real tCoefficients[ 6 ][ 3 ] = {};

        for ( uint iVertex = 0; iVertex < 8; ++iVertex )
        {
            const real tXi           = tSigns[ iVertex ][ 0 ];
            const real tEta          = tSigns[ iVertex ][ 1 ];
            const real tZeta         = tSigns[ iVertex ][ 2 ];
            const real tWeights[ 6 ] = { tXi, tEta, tZeta, tXi * tEta, tEta * tZeta, tXi * tZeta };

            for ( uint iCoefficient = 0; iCoefficient < 6; ++iCoefficient )
            {
                for ( uint iDirection = 0; iDirection < 3; ++iDirection )
                {
                    tCoefficients[ iCoefficient ][ iDirection ] += 0.125 * tWeights[ iCoefficient ]
                                                                 * ( tCoordinates( iVertex, iDirection ) - tCoordinates( 0, iDirection ) );
                }
            }
        }

        const uint tFirstDirection  = ( aDirection + 1 ) % 3;
        const uint tSecondDirection = ( aDirection + 2 ) % 3;
        auto       tCrossComponent  = [ & ]( uint aFirst, uint aSecond ) {
            return tCoefficients[ aFirst ][ tFirstDirection ] * tCoefficients[ aSecond ][ tSecondDirection ]
                 - tCoefficients[ aFirst ][ tSecondDirection ] * tCoefficients[ aSecond ][ tFirstDirection ];
        };

        // Integrating det(J) exactly gives
        // V = 8*det(a1,a2,a3) + 8/3*(det(a1,a4,a6) + det(a4,a2,a5) + det(a6,a5,a3)).
        // Differentiate this polynomial using da_k/dx_vertex = sign_product/8.
        // Array indices 0,...,5 represent a1,...,a6; cyclic components give each cross product.
        const real tXi   = tSigns[ aLocalVertexID ][ 0 ];
        const real tEta  = tSigns[ aLocalVertexID ][ 1 ];
        const real tZeta = tSigns[ aLocalVertexID ][ 2 ];

        return tXi * ( tCrossComponent( 1, 2 ) + tCrossComponent( 3, 5 ) / 3.0 )
             + tEta * ( tCrossComponent( 2, 0 ) + tCrossComponent( 4, 3 ) / 3.0 )
             + tZeta * ( tCrossComponent( 0, 1 ) + tCrossComponent( 5, 4 ) / 3.0 )
             + tXi * tEta * ( tCrossComponent( 5, 0 ) + tCrossComponent( 1, 4 ) ) / 3.0
             + tEta * tZeta * ( tCrossComponent( 3, 1 ) + tCrossComponent( 2, 5 ) ) / 3.0
             + tXi * tZeta * ( tCrossComponent( 0, 3 ) + tCrossComponent( 4, 2 ) ) / 3.0;
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
