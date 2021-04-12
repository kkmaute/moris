/*
 * cl_MTK_Cell_Info_Hex27.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Hex27.hpp"
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
        Cell_Info_Hex27::get_cell_geometry() const
        {
            return Geometry_Type::HEX;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Hex27::get_cell_interpolation_order() const
        {
            return Interpolation_Order::QUADRATIC;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Hex27::get_cell_integration_order() const
        {
            return Integration_Order::HEX_4x4x4;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Hex27::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // init cell shape
            CellShape tCellShape = CellShape::RECTANGULAR;

            // error threshold
            real tEpsilon = 1.0E-8;

            // looping through each face
            for ( uint iFace = 0; iFace < 6; iFace++)
            {
                // getting nodes on the face
                moris::Matrix<moris::IndexMat> tFaceNodes = this->get_node_to_face_map( iFace );

                // get edges to define check plane
                moris::Matrix< DDRMat > tEdge0 = tVertices( tFaceNodes( 1 ) )->get_coords() - tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge1 = tVertices( tFaceNodes( 2 ) )->get_coords() - tVertices( tFaceNodes( 1 ) )->get_coords();
                moris::Matrix< DDRMat > tFaceVec01 = cross(tEdge0,tEdge1);

                // these are the other plane normals produced by the other 6 nodes.
                // keeping node 0 to ensure the plane is at the same offset.
                moris::Matrix< DDRMat > tEdge2 = tVertices( tFaceNodes( 3 ) )->get_coords() - tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge3 = tVertices( tFaceNodes( 4 ) )->get_coords() - tVertices( tFaceNodes( 3 ) )->get_coords();
                moris::Matrix< DDRMat > tFaceVec23 = cross(tEdge2,tEdge3);

                moris::Matrix< DDRMat > tEdge4 = tVertices( tFaceNodes( 5 ) )->get_coords() - tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge5 = tVertices( tFaceNodes( 6 ) )->get_coords() - tVertices( tFaceNodes( 5 ) )->get_coords();
                moris::Matrix< DDRMat > tFaceVec45 = cross(tEdge4,tEdge5);

                moris::Matrix< DDRMat > tEdge6 = tVertices( tFaceNodes( 7 ) )->get_coords() - tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge7 = tVertices( tFaceNodes( 8 ) )->get_coords() - tVertices( tFaceNodes( 7 ) )->get_coords();
                moris::Matrix< DDRMat > tFaceVec67 = cross(tEdge6,tEdge7);

                // All three of the plane normals must be parallel in order to be considered straight shape
                if( norm( cross(tFaceVec01, tFaceVec23) ) > tEpsilon ||
                    norm( cross(tFaceVec01, tFaceVec45) ) > tEpsilon ||
                    norm( cross(tFaceVec01, tFaceVec67) ) > tEpsilon )
                {
                    tCellShape = CellShape::GENERAL;
                    break;
                }

                // Checking if this face is perpindicular to the next face plane
                
                if( iFace < 5 && tCellShape == CellShape::RECTANGULAR )
                {
                    moris::Matrix<moris::IndexMat> tNextFaceNodes = this->get_node_to_face_map( iFace + 1 );

                    // get the plane defined by the first 3 nodes on this face
                    moris::Matrix< DDRMat > tNextEdge0 = tVertices( tNextFaceNodes( 1 ) )->get_coords() -
                            tVertices( tNextFaceNodes( 0 ) )->get_coords();
                    moris::Matrix< DDRMat > tNextEdge1 = tVertices( tNextFaceNodes( 2 ) )->get_coords() -
                            tVertices( tNextFaceNodes( 1 ) )->get_coords();
                    moris::Matrix< DDRMat > tNextFaceVec01 = cross(tNextEdge0,tNextEdge1);

                    if( iFace < 4 )
                    {
                        // the 2 faces should be perpindicular in order for the shape to be rectangular
                        if( std::abs( dot(tFaceVec01, tNextFaceVec01) ) > tEpsilon )
                        {
                            tCellShape = CellShape::STRAIGHT;
                        }
                    }
                    // the last 2 faces should be parallel in order for it to be 
                    else
                    {
                        // the 2 faces should be parallel in order for the shape to be rectangular
                        if( norm( cross(tFaceVec01, tNextFaceVec01) ) > tEpsilon )
                        {
                            tCellShape = CellShape::STRAIGHT;
                        }
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex27::get_num_verts() const
        {
            return 27;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex27::get_num_facets() const
        {
            return 6;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex27::get_num_edges() const
        {
            return 12;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex27::get_num_verts_per_facet() const
        {
            return 9;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex27::get_loc_coord_dim() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        inline
        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_face_map() const
        {
            return {{0, 1, 5, 4,  8, 13, 16, 12, 25},
                {1, 2, 6, 5,  9, 14, 17, 13, 24},
                {2, 3, 7, 6, 10, 15, 18, 14, 26},
                {0, 4, 7, 3, 12, 19, 15, 11, 23},
                {0, 3, 2, 1, 11, 10,  9,  8, 21},
                {4, 5, 6, 7, 16, 17, 18, 19, 22}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_edge_map() const
        {
            return {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};
        }

        // ----------------------------------------------------------------------------------

        inline
        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{0, 1, 5, 4,  8, 13, 16, 12, 25}}; break; }
                case(1):{ return {{1, 2, 6, 5,  9, 14, 17, 13, 24}}; break; }
                case(2):{ return {{2, 3, 7, 6, 10, 15, 18, 14, 26}}; break; }
                case(3):{ return {{0, 4, 7, 3, 12, 19, 15, 11, 23}}; break; }
                case(4):{ return {{0, 3, 2, 1, 11, 10,  9,  8, 21}}; break; }
                case(5):{ return {{4, 5, 6, 7, 16, 17, 18, 19, 22}}; break; }
                default:
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
            }
        }

        // ----------------------------------------------------------------------------------
        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
                case( 0):{ return {{0, 1}}; break; }
                case( 1):{ return {{1, 2}}; break; }
                case( 2):{ return {{2, 3}}; break; }
                case( 3):{ return {{3, 0}}; break; }
                case( 4):{ return {{4, 5}}; break; }
                case( 5):{ return {{5, 6}}; break; }
                case( 6):{ return {{6, 7}}; break; }
                case( 7):{ return {{7, 4}}; break; }
                case( 8):{ return {{0, 4}}; break; }
                case( 9):{ return {{1, 5}}; break; }
                case(10):{ return {{2, 6}}; break; }
                case(11):{ return {{3, 7}}; break; }
                default:
                    MORIS_ASSERT(0,"Invalid edge ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_geometric_node_to_facet_map() const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::uint
        Cell_Info_Hex27::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case 0:{ return 2; break; }
                case 1:{ return 3; break; }
                case 2:{ return 0; break; }
                case 3:{ return 1; break; }
                case 4:{ return 5; break; }
                case 5:{ return 4; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return MORIS_UINT_MAX;
                    break;
                }
            }
        }

        Matrix<DDRMat>
        Cell_Info_Hex27::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
                case  0: { return {{ -1.0, -1.0, -1.0 }}; break; }
                case  1: { return {{ +1.0, -1.0, -1.0 }}; break; }
                case  2: { return {{ +1.0, +1.0, -1.0 }}; break; }
                case  3: { return {{ -1.0, +1.0, -1.0 }}; break; }
                case  4: { return {{ -1.0, -1.0, +1.0 }}; break; }
                case  5: { return {{ +1.0, -1.0, +1.0 }}; break; }
                case  6: { return {{ +1.0, +1.0, +1.0 }}; break; }
                case  7: { return {{ -1.0, +1.0, +1.0 }}; break; }
                case  8: { return {{ +0.0, -1.0, -1.0 }}; break; }
                case  9: { return {{ +1.0,  0.0, -1.0 }}; break; }
                case 10: { return {{  0.0, +1.0, -1.0 }}; break; }
                case 11: { return {{ -1.0,  0.0, -1.0 }}; break; }
                case 12: { return {{ -1.0, -1.0,  0.0 }}; break; }
                case 13: { return {{ +1.0, -1.0,  0.0 }}; break; }
                case 14: { return {{ +1.0, +1.0,  0.0 }}; break; }
                case 15: { return {{ -1.0, +1.0,  0.0 }}; break; }
                case 16: { return {{  0.0, -1.0, +1.0 }}; break; }
                case 17: { return {{ +1.0,  0.0, +1.0 }}; break; }
                case 18: { return {{  0.0, +1.0, +1.0 }}; break; }
                case 19: { return {{ -1.0,  0.0, +1.0 }}; break; }
                case 20: { return {{  0.0,  0.0,  0.0 }}; break; }
                case 21: { return {{  0.0,  0.0, -1.0 }}; break; }
                case 22: { return {{  0.0,  0.0, +1.0 }}; break; }
                case 23: { return {{ -1.0,  0.0,  0.0 }}; break; }
                case 24: { return {{ +1.0,  0.0,  0.0 }}; break; }
                case 25: { return {{  0.0, -1.0,  0.0 }}; break; }
                case 26: { return {{  0.0, +1.0,  0.0 }}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid vertex ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------
        inline
        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex27::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{0,1},{1,5}}; break; }
                case(1):{ return {{1,2},{2,6}}; break; }
                case(2):{ return {{2,3},{3,7}}; break; }
                case(3):{ return {{7,3},{3,0}}; break; }
                case(4):{ return {{0,3},{3,2}}; break; }
                case(5):{ return {{4,5},{5,6}}; break; }
                default:
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real
         Cell_Info_Hex27::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // FIXME: only works for rectangular cells
            Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
            Matrix<DDRMat> tNode6Coords = tVertices(6)->get_coords();

            real tLx = std::abs(tNode0Coords(0) - tNode6Coords(0));
            real tLy = std::abs(tNode0Coords(1) - tNode6Coords(1));
            real tLz = std::abs(tNode0Coords(2) - tNode6Coords(2));

            return tLx*tLy*tLz;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex27::compute_cell_side_size( moris::mtk::Cell const * aCell ,
                moris_index const & aSideOrd) const
        {
            moris::Cell< Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            // FIXME: only works for rectangular cells
            Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords1 = tVertices(1)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(3)->get_coords();

            return norm( cross( tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0 ) );
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Hex27::eval_N( const Matrix< DDRMat > & aXi,
                Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX27 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            auto    xi = aXi( 0 );
            auto   eta = aXi( 1 );
            auto  zeta = aXi( 2 );

            // often used constants
            auto    xi2 = std::pow(   xi, 2 );
            auto   eta2 = std::pow(  eta, 2 );
            auto  zeta2 = std::pow( zeta, 2 );

            auto a = -0.25 * eta * zeta;
            auto b = -0.25 * xi * zeta;
            auto c = -0.25 * xi * eta;
            auto d = 0.125 * xi * eta * zeta;

            // populate output matrix
            aNXi.set_size(1,27);
            aNXi(  0 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aNXi(  1 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aNXi(  2 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aNXi(  3 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aNXi(  4 ) = d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aNXi(  5 ) = d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aNXi(  6 ) = d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aNXi(  7 ) = d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aNXi(  8 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
            aNXi(  9 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
            aNXi( 10 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
            aNXi( 11 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
            aNXi( 12 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
            aNXi( 13 ) = c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
            aNXi( 14 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
            aNXi( 15 ) = c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
            aNXi( 16 ) = a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
            aNXi( 17 ) = b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
            aNXi( 18 ) = a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
            aNXi( 19 ) = b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
            aNXi( 20 ) = -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
            aNXi( 21 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
            aNXi( 22 ) = ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
            aNXi( 23 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
            aNXi( 24 ) = ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
            aNXi( 25 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
            aNXi( 26 ) = ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;
        }
    }
}
