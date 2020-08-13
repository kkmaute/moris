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
        Cell_Info_Quad9::compute_cell_size( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            Matrix<DDRMat> tNode1Coords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();

            real tLx = std::abs(tNode1Coords0(0) - tNodeCoords2(0));
            real tLy = std::abs(tNode1Coords0(1) - tNodeCoords2(1));

            return tLx*tLy;
        }
        // ----------------------------------------------------------------------------------
        moris::real
        Cell_Info_Quad9::compute_cell_side_size( moris::mtk::Cell const * aCell ,
                moris_index const & aSideOrd) const
        {
            moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            Matrix<DDRMat> tLVec = tVertices(1)->get_coords() - tVertices(0)->get_coords();

            return moris::norm(tLVec);
        }
        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Quad9::eval_N( const Matrix< DDRMat > & aXi,
                                       Matrix< DDRMat > & aNXi ) const
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


