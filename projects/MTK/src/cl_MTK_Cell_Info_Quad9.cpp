/*
 * cl_MTK_Cell_Info_Quad9.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Quad9.hpp"
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
#include "fn_cross.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------

        enum Geometry_Type
        Cell_Info_Quad9::get_cell_geometry() const
        {
            return Geometry_Type::QUAD;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Quad9::get_cell_interpolation_order() const
        {
            return Interpolation_Order::QUADRATIC;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Quad9::get_cell_integration_order() const
        {
            return Integration_Order::QUAD_3x3;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Quad9::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // init cell shape
            CellShape tCellShape = CellShape::RECTANGULAR;

            // error threshold
            real tEpsilon = 1.0E-8;

            // looping through each Edge
            for ( uint iEdge = 0; iEdge < 4; iEdge++)
            {
                // getting nodes on the edge
                moris::Matrix<moris::IndexMat> tEdgeNodes = this->get_node_to_edge_map( iEdge );

                // getting getting vectors on this edge
                moris::Matrix< DDRMat > tEdgeVec0 = tVertices( tEdgeNodes(1) )->get_coords() - tVertices( tEdgeNodes(0) )->get_coords();
                moris::Matrix< DDRMat > tEdgeVec1 = tVertices( tEdgeNodes(2) )->get_coords() - tVertices( tEdgeNodes(1) )->get_coords();

                // perform cross product of the 2D vectors
                real tCross = tEdgeVec0(0) * tEdgeVec1(1) - tEdgeVec0(1) * tEdgeVec1(0);

                // are the nodes on this edge on the same vector?
                if ( std::abs( tCross ) > tEpsilon )
                {
                    tCellShape = CellShape::GENERAL;
                    break;
                }
 
                // checking if this edge is perpindicular to the next edge ordinal. Only need to do it for first 3 corners
                if ( iEdge < 3 )
                {
                    // next adjacent edge. Using the mid edge node for convenience
                    moris::Matrix< DDRMat > tEdgeVec2 = tVertices( iEdge + 5 )->get_coords() - tVertices( iEdge + 1 )->get_coords();
                    
                    if ( std::abs( dot( tEdgeVec0, tEdgeVec2 ) ) > tEpsilon )
                    {
                        tCellShape = CellShape::STRAIGHT;
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad9::get_num_verts() const
        {
            return 9;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad9::get_num_facets() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad9::get_num_edges() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad9::get_num_verts_per_facet() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad9::get_loc_coord_dim() const
        {
            return 2;
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_to_face_map() const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_to_edge_map() const
        {
            return {{0,1,4}, {1,2,5}, {2,3,6}, {3,0,7}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_to_facet_map() const
        {
            return this->get_node_to_edge_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_geometric_node_to_facet_map() const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.get_node_to_facet_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            Cell_Info_Quad4 tQuad4;
            return tQuad4.get_node_to_facet_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
                case(0):{ return {{0,1,4}}; break; }
                case(1):{ return {{1,2,5}}; break; }
                case(2):{ return {{2,3,6}}; break; }
                case(3):{ return {{3,0,7}}; break; }
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
        Cell_Info_Quad9::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad9::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{1,0}}; break; }
                case(1):{ return {{2,1}}; break; }
                case(2):{ return {{3,2}}; break; }
                case(3):{ return {{0,3}}; break; }
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
        Cell_Info_Quad9::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
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

        Matrix<DDRMat>
        Cell_Info_Quad9::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
                case 0: { return {{-1.000000000000000e+00,  -1.000000000000000e+00}}; break; }
                case 1: { return {{+1.000000000000000e+00,  -1.000000000000000e+00}}; break; }
                case 2: { return {{+1.000000000000000e+00,  +1.000000000000000e+00}}; break; }
                case 3: { return {{-1.000000000000000e+00,  +1.000000000000000e+00}}; break; }
                case 4: { return {{ 0.000000000000000e+00,  -1.000000000000000e+00}}; break; }
                case 5: { return {{+1.000000000000000e+00,   0.000000000000000e+00}}; break; }
                case 6: { return {{ 0.000000000000000e+00,  +1.000000000000000e+00}}; break; }
                case 7: { return {{-1.000000000000000e+00,   0.000000000000000e+00}}; break; }
                case 8: { return {{ 0.000000000000000e+00,   0.000000000000000e+00}}; break; }
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
        Cell_Info_Quad9::get_loc_coord_on_side_ordinal(moris::uint aSideOrdinal) const
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
        Cell_Info_Quad9::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2,"Cell_Info_Quad4::compute_cell_size_special only works in 2D.\n");

            // FIXME: only works for rectangular cells
            real tLx = std::abs(tNodeCoords0(0) - tNodeCoords2(0));
            real tLy = std::abs(tNodeCoords0(1) - tNodeCoords2(1));

            return tLx*tLy;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad9::compute_cell_side_size(
                moris::mtk::Cell const * aCell ,
                moris_index      const & aSideOrd) const
        {
            moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            Matrix<DDRMat> tLVec = tVertices(1)->get_coords() - tVertices(0)->get_coords();

            return moris::norm(tLVec);
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Quad9::eval_N(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2,
                    "QUAD9 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // often used constants
            real    c = xi * eta * 0.25;
            real  xi2 = std::pow(  xi, 2 );
            real eta2 = std::pow( eta, 2 );

            // populate output matrix
            aNXi.set_size( 1, 9 );
            aNXi( 0 ) = ( c * ( eta - 1.0 ) * (xi - 1.0) );
            aNXi( 1 ) = ( c * ( eta - 1.0 ) * (xi + 1.0) );
            aNXi( 2 ) = ( c * ( eta + 1.0 ) * (xi + 1.0) );
            aNXi( 3 ) = ( c * ( eta + 1.0 ) * (xi - 1.0) );
            aNXi( 4 ) = ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
            aNXi( 5 ) = ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
            aNXi( 6 ) = ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
            aNXi( 7 ) = ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
            aNXi( 8 ) = ( eta2 - 1.0 )*( xi2 - 1.0 );
        }
    }
}
