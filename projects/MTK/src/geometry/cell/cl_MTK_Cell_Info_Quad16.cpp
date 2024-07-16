/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Quad16.cpp
 *
 */

#include "cl_MTK_Cell_Info_Quad16.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
#include "op_times.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------

        enum Geometry_Type
        Cell_Info_Quad16::get_cell_geometry() const
        {
            return Geometry_Type::QUAD;
        }

        // ----------------------------------------------------------------------------------

        enum CellTopology
        Cell_Info_Quad16::get_cell_topology() const
        {
            return CellTopology::QUAD16;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Quad16::get_cell_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Quad16::get_cell_integration_order() const
        {
            return Integration_Order::QUAD_5x5;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Quad16::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            // throw error, since this function hasn't been tested. Would need to be able to create a cell without an STK mesh
//            MORIS_ERROR(false,"Cell_Info_Quad16::compute_cell_shape - this has not been tested. no Quad16 STK implementation");

            // getting vertices and storing them in a local matrix, since each node will be used a few times
            Vector< Vertex* > tVertices = aCell->get_vertex_pointers();

            // init cell shape
            CellShape tCellShape = CellShape::RECTANGULAR;

            // error threshold
            real tEpsilon = 1.0E-8;

            // init cell of edge vectors
            Vector< moris::Matrix< DDRMat > > tEdgeVectors( 4 );

            // looping through each Edge
            for ( uint iEdge = 0; iEdge < 4; iEdge++)
            {
                // getting nodes on the edge
                moris::Matrix<moris::IndexMat> tEdgeNodes = this->get_node_to_edge_map( iEdge );

                // get node 0 coordinates
                moris::Matrix< DDRMat > tVertex0 = tVertices( tEdgeNodes(0) )->get_coords();

                // getting getting vectors on this edge
                tEdgeVectors( iEdge )             = tVertices( tEdgeNodes(1) )->get_coords() - tVertex0;
                moris::Matrix< DDRMat > tEdgeVec1 = tVertices( tEdgeNodes(2) )->get_coords() - tVertex0;
                moris::Matrix< DDRMat > tEdgeVec2 = tVertices( tEdgeNodes(3) )->get_coords() - tVertex0;

                // perform cross product of the 2D vectors
                auto tCross01 = tEdgeVectors( iEdge )(0) * tEdgeVec1(1) - tEdgeVectors( iEdge )(1) * tEdgeVec1(0);
                auto tCross12 = tEdgeVec1(0) * tEdgeVec2(1)             - tEdgeVec1(1) * tEdgeVec2(0);

                // are the nodes on this edge on the same vector?
                if ( std::abs( tCross01 ) > tEpsilon ||
                     std::abs( tCross12 ) > tEpsilon  )
                {
                    tCellShape = CellShape::GENERAL;
                    break;
                }
            }

            // check if the cell is still rectangular check if it is parallel or straight
            if ( tCellShape == CellShape::RECTANGULAR )
            {
                // get cross products of opposite edges
                auto tCross02 = tEdgeVectors(0)(0)*tEdgeVectors(2)(1) - tEdgeVectors(0)(1)*tEdgeVectors(2)(0);
                auto tCross13 = tEdgeVectors(1)(0)*tEdgeVectors(3)(1) - tEdgeVectors(1)(1)*tEdgeVectors(3)(0);

                // if they aren't parallel, then it is a straight cell
                if( std::abs( tCross02 ) > tEpsilon ||
                    std::abs( tCross13 ) > tEpsilon )
                {
                    tCellShape = CellShape::STRAIGHT;
                }

                // if they are parallel edges check if it is rectangular and aligned
                else
                {
                    // if the first two edges aren't perpendicular or if it isn't aligned
                    if( std::abs(dot( tEdgeVectors(0), tEdgeVectors(1) ) ) > tEpsilon ||
                        std::abs( tEdgeVectors(0)(1) )                     > tEpsilon ||
                        tEdgeVectors(0)(0)                                 < 0.0 )
                    {
                        tCellShape = CellShape::PARALLEL;
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad16::get_num_verts() const
        {
            return 16;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad16::get_num_facets() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad16::get_num_edges() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad16::get_num_verts_per_facet() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad16::get_loc_coord_dim() const
        {
            return 2;
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_face_map() const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_edge_map() const
        {
            return {
                { 0, 1,  4,  5 },
                { 1, 2,  6,  7 },
                { 2, 3,  8,  9 },
                { 3, 0, 10, 11 } };
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_facet_map() const
        {
            return this->get_node_to_edge_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_geometric_node_to_facet_map() const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.get_node_to_facet_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.get_node_to_facet_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------
        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
                case(0):{ return {{ 0, 1,  4,  5 }}; break; }
                case(1):{ return {{ 1, 2,  6,  7 }}; break; }
                case(2):{ return {{ 2, 3,  8,  9 }}; break; }
                case(3):{ return {{ 3, 0, 10, 11 }}; break; }
                default:
                {
                    MORIS_ASSERT(0,"Invalid edge ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad16::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{ 1, 0 }}; break; }
                case(1):{ return {{ 2, 1 }}; break; }
                case(2):{ return {{ 3, 2 }}; break; }
                case(3):{ return {{ 0, 3 }}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::uint
        Cell_Info_Quad16::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return 2; break; }
                case(1):{ return 3; break; }
                case(2):{ return 0; break; }
                case(3):{ return 1; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return MORIS_UINT_MAX;
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        Vector< moris_index >
        Cell_Info_Quad16::get_vertex_path_to_entity_rank_and_ordinal(
            moris_index aVertexOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
        {

            switch ( aOtherEntityRank )
            {
                // node to node paths
                case 0:
                {
                    const Matrix< IndexMat > tVertexToVertexRanks = {
                        { -1, +1, +3, +1, +1, +1, +3, +3, +3, +3, +1, +1, +3, +3, +3, +3},
                        { +1, -1, +1, +3, +1, +1, +1, +1, +3, +3, +3, +3, +3, +3, +3, +3},
                        { +3, +1, -1, +1, +3, +3, +1, +1, +1, +1, +3, +3, +3, +3, +3, +3},
                        { +1, +3, +1, -1, +3, +3, +3, +3, +1, +1, +1, +1, +3, +3, +3, +3},
                        { +1, +1, +3, +3, -1, +1, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3},
                        { +1, +1, +3, +3, +1, -1, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3},
                        { +3, +1, +1, +3, +3, +3, -1, +1, +3, +3, +3, +3, +3, +3, +3, +3},
                        { +3, +1, +1, +3, +3, +3, +1, -1, +3, +3, +3, +3, +3, +3, +3, +3},
                        { +3, +3, +1, +1, +3, +3, +3, +3, -1, +1, +3, +3, +3, +3, +3, +3},
                        { +3, +3, +1, +1, +3, +3, +3, +3, +1, -1, +3, +3, +3, +3, +3, +3},
                        { +1, +3, +3, +1, +3, +3, +3, +3, +3, +3, -1, +1, +3, +3, +3, +3},
                        { +1, +3, +3, +1, +3, +3, +3, +3, +3, +3, +1, -1, +3, +3, +3, +3},
                        { +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, -1, +3, +3, +3},
                        { +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, -1, +3, +3},
                        { +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, -1, +3},
                        { +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, +3, -1} };

                    const Matrix< IndexMat > tVertexToVertexIndices = {
                        { -1, +0, +0, +3, +0, +0, +0, +0, +0, +0, +3, +3, +0, +0, +0, +0},
                        { +0, -1, +1, +0, +0, +0, +1, +1, +0, +0, +0, +0, +0, +0, +0, +0},
                        { +0, +1, -1, +2, +0, +0, +1, +1, +2, +2, +0, +0, +0, +0, +0, +0},
                        { +3, +0, +2, -1, +0, +0, +0, +0, +2, +2, +3, +3, +0, +0, +0, +0},
                        { +0, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0},
                        { +0, +0, +0, +0, +0, -1, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0},
                        { +0, +1, +1, +0, +0, +0, -1, +1, +0, +0, +0, +0, +0, +0, +0, +0},
                        { +0, +1, +1, +0, +0, +0, +1, -1, +0, +0, +0, +0, +0, +0, +0, +0},
                        { +0, +0, +2, +2, +0, +0, +0, +0, -1, +2, +0, +0, +0, +0, +0, +0},
                        { +0, +0, +2, +2, +0, +0, +0, +0, +2, -1, +0, +0, +0, +0, +0, +0},
                        { +3, +0, +0, +3, +0, +0, +0, +0, +0, +0, -1, +3, +0, +0, +0, +0},
                        { +3, +0, +0, +3, +0, +0, +0, +0, +0, +0, +3, -1, +0, +0, +0, +0},
                        { +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, -1, +0, +0, +0},
                        { +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, -1, +0, +0},
                        { +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, -1, +0},
                        { +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, +0, -1} };

                    moris_index tPathRank = tVertexToVertexRanks( (uint) aVertexOrdinal, (uint) aOtherEntityOrdinal );
                    moris_index tPathIndex = tVertexToVertexIndices( (uint) aVertexOrdinal, (uint) aOtherEntityOrdinal );

                    MORIS_ASSERT( tPathRank != -1 && tPathIndex != -1,
                        "Cell_Info_Quad16::get_vertex_path_to_entity_rank_and_ordinal() - Vertex doesn't have path to itself." );

                    return { tPathIndex, tPathRank };
                    break;
                }

                // node to edge paths
                case 1:
                {
                    const Matrix< IndexMat > tVertexToEdgeRanks = {
                        { +1, +3, +3, +1},
                        { +1, +1, +3, +3},
                        { +3, +1, +1, +3},
                        { +3, +3, +1, +1},
                        { +1, +3, +3, +3},
                        { +1, +3, +3, +3},
                        { +3, +1, +3, +3},
                        { +3, +1, +3, +3},
                        { +3, +3, +1, +3},
                        { +3, +3, +1, +3},
                        { +3, +3, +3, +1},
                        { +3, +3, +3, +1},
                        { +3, +3, +3, +3},
                        { +3, +3, +3, +3},
                        { +3, +3, +3, +3},
                        { +3, +3, +3, +3} };

                    const Matrix< IndexMat > tVertexToEdgeIndices = {
                        { +0, +0, +0, +3},
                        { +0, +1, +0, +0},
                        { +0, +1, +2, +0},
                        { +0, +0, +2, +3},
                        { +0, +0, +0, +0},
                        { +0, +0, +0, +0},
                        { +0, +1, +0, +0},
                        { +0, +1, +0, +0},
                        { +0, +0, +2, +0},
                        { +0, +0, +2, +0},
                        { +0, +0, +0, +3},
                        { +0, +0, +0, +3},
                        { +0, +0, +0, +0},
                        { +0, +0, +0, +0},
                        { +0, +0, +0, +0},
                        { +0, +0, +0, +0} };

                    moris_index tPathRank = tVertexToEdgeRanks( (uint) aVertexOrdinal, (uint) aOtherEntityOrdinal );
                    moris_index tPathIndex = tVertexToEdgeIndices( (uint) aVertexOrdinal, (uint) aOtherEntityOrdinal );

                    MORIS_ASSERT( tPathRank != -1 && tPathIndex != -1,
                        "Cell_Info_Quad16::get_vertex_path_to_entity_rank_and_ordinal() - Vertex doesn't have path to itself." );

                    return { tPathIndex, tPathRank };
                    break;
                }

                default:
                {
                    MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                    return Vector< moris_index >( 0 );
                }
            } // end: switch aOtherEntityRank
        }

        // ----------------------------------------------------------------------------------

        Vector< moris_index >
        Cell_Info_Quad16::get_edge_path_to_entity_rank_and_ordinal(
            moris_index aEdgeOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.get_edge_path_to_entity_rank_and_ordinal( aEdgeOrdinal, aOtherEntityOrdinal, aOtherEntityRank );
        }

        // ----------------------------------------------------------------------------------

        bool
        Cell_Info_Quad16::is_entity_connected_to_facet(
            moris_index aFacetOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.is_entity_connected_to_facet( aFacetOrdinal, aOtherEntityOrdinal, aOtherEntityRank );
        }

        // ----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Cell_Info_Quad16::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
                case  0:{ return {{ -1.000000000000000e+00,  -1.000000000000000e+00 }}; break;}
                case  1:{ return {{ +1.000000000000000e+00,  -1.000000000000000e+00 }}; break;}
                case  2:{ return {{ +1.000000000000000e+00,  +1.000000000000000e+00 }}; break;}
                case  3:{ return {{ -1.000000000000000e+00,  +1.000000000000000e+00 }}; break;}
                case  4:{ return {{ -3.333333333333334e-01,  -1.000000000000000e+00 }}; break;}
                case  5:{ return {{ +3.333333333333333e-01,  -1.000000000000000e+00 }}; break;}
                case  6:{ return {{ +1.000000000000000e+00,  -3.333333333333334e-01 }}; break;}
                case  7:{ return {{ +1.000000000000000e+00,  +3.333333333333333e-01 }}; break;}
                case  8:{ return {{ +3.333333333333333e-01,  +1.000000000000000e+00 }}; break;}
                case  9:{ return {{ -3.333333333333334e-01,  +1.000000000000000e+00 }}; break;}
                case 10:{ return {{ -1.000000000000000e+00,  +3.333333333333333e-01 }}; break;}
                case 11:{ return {{ -1.000000000000000e+00,  -3.333333333333334e-01 }}; break;}
                case 12:{ return {{ -3.333333333333334e-01,  -3.333333333333334e-01 }}; break;}
                case 13:{ return {{ +3.333333333333333e-01,  -3.333333333333334e-01 }}; break;}
                case 14:{ return {{ +3.333333333333333e-01,  +3.333333333333333e-01 }}; break;}
                case 15:{ return {{ -3.333333333333334e-01,  +3.333333333333333e-01 }}; break;}
                default:
                {
                    MORIS_ERROR(0,"Invalid vertex ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------
        Matrix<DDRMat>
        Cell_Info_Quad16::get_loc_coord_on_side_ordinal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{-1,-1 }, { 1,-1 }}; break; }
                case(1):{ return {{ 1,-1 }, { 1, 1 }}; break; }
                case(2):{ return {{ 1, 1 }, {-1, 1 }}; break; }
                case(3):{ return {{-1, 1 }, {-1,-1 }}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad16::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            Vector< Vertex* > tVertices = aCell->get_vertex_pointers();

            Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2,"Cell_Info_Quad16::compute_cell_size_special only works in 2D.\n");

            // FIXME: only works for rectangular cells
            real tLx = std::abs(tNodeCoords0(0) - tNodeCoords2(0));
            real tLy = std::abs(tNodeCoords0(1) - tNodeCoords2(1));

            return tLx*tLy;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad16::compute_cell_side_size( moris::mtk::Cell const * aCell ,
                moris_index const & aSideOrd) const
        {
            Vector< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            Matrix<DDRMat> tLVec = tVertices(1)->get_coords() - tVertices(0)->get_coords();

            return moris::norm(tLVec);
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Quad16::eval_N( const Matrix< DDRMat > & aXi,
                Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD16 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

            // populate matrix with values
            aNXi.set_size(1,16);
            aNXi(  0 ) = a0*b0;
            aNXi(  1 ) = a3*b0;
            aNXi(  2 ) = a3*b3;
            aNXi(  3 ) = a0*b3;
            aNXi(  4 ) = a1*b0;
            aNXi(  5 ) = a2*b0;
            aNXi(  6 ) = a3*b1;
            aNXi(  7 ) = a3*b2;
            aNXi(  8 ) = a2*b3;
            aNXi(  9 ) = a1*b3;
            aNXi( 10 ) = a0*b2;
            aNXi( 11 ) = a0*b1;
            aNXi( 12 ) = a1*b1;
            aNXi( 13 ) = a2*b1;
            aNXi( 14 ) = a2*b2;
            aNXi( 15 ) = a1*b2;
        }
    }
}

