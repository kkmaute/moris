/*
 * cl_MTK_Cell_Info_Hex8.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

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
        Cell_Info_Hex8::get_cell_geometry() const
        {
            return Geometry_Type::HEX;
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
        Cell_Info_Hex8::compute_cell_shape(moris::mtk::Cell const *aCell) const
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

                moris::Matrix< DDRMat > tEdge0 = tVertices( tFaceNodes( 1 ) )->get_coords() - tVertices( tFaceNodes( 0 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge1 = tVertices( tFaceNodes( 2 ) )->get_coords() - tVertices( tFaceNodes( 1 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge2 = tVertices( tFaceNodes( 3 ) )->get_coords() - tVertices( tFaceNodes( 2 ) )->get_coords();
                moris::Matrix< DDRMat > tEdge3 = tVertices( tFaceNodes( 0 ) )->get_coords() - tVertices( tFaceNodes( 3 ) )->get_coords();

                // are the edges perpindicular to the adjacent edge?
                if ( std::abs( dot(tEdge1, tEdge0 ) ) > tEpsilon || 
                     std::abs( dot(tEdge2, tEdge1 ) ) > tEpsilon ||
                     std::abs( dot(tEdge3, tEdge2 ) ) > tEpsilon ||
                     std::abs( dot(tEdge0, tEdge3 ) ) > tEpsilon )
                {
                    // cross products at opposite corners
                    moris::Matrix< DDRMat > tCross0 = cross( tEdge1, tEdge0 );
                    moris::Matrix< DDRMat > tCross1 = cross( tEdge3, tEdge2 );

                    // if the the cross products at the corners are not parallel
                    if ( norm( cross( tCross0, tCross1 )) > tEpsilon )
                    {
                        // no need to check other faces if this is true
                        tCellShape = CellShape::GENERAL;
                        break;
                    }

                    // setting cell shape to straight if the face is planar
                    else
                    {
                        tCellShape = CellShape::STRAIGHT;
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

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_face_map() const
        {
            return {{1, 5, 4, 0}, {1, 2, 6, 5}, {3, 7, 6, 2}, {0, 4, 7, 3}, {0, 3, 2, 1}, {4, 5, 6, 7}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_edge_map() const
        {
            return {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
            case (0):
            {
                return {{1, 5, 4, 0}};
                break;
            }
            case (1):
            {
                return {{1, 2, 6, 5}};
                break;
            }
            case (2):
            {
                return {{3, 7, 6, 2}};
                break;
            }
            case (3):
            {
                return {{0, 4, 7, 3}};
                break;
            }
            case (4):
            {
                return {{0, 3, 2, 1}};
                break;
            }
            case (5):
            {
                return {{4, 5, 6, 7}};
                break;
            }
            default:
                MORIS_ERROR(0, "Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0, 0);
                break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
            case (0):
            {
                return {{0, 1}};
                break;
            }
            case (1):
            {
                return {{1, 2}};
                break;
            }
            case (2):
            {
                return {{2, 3}};
                break;
            }
            case (3):
            {
                return {{3, 0}};
                break;
            }
            case (4):
            {
                return {{4, 5}};
                break;
            }
            case (5):
            {
                return {{5, 6}};
                break;
            }
            case (6):
            {
                return {{6, 7}};
                break;
            }
            case (7):
            {
                return {{7, 4}};
                break;
            }
            case (8):
            {
                return {{0, 4}};
                break;
            }
            case (9):
            {
                return {{1, 5}};
                break;
            }
            case (10):
            {
                return {{2, 6}};
                break;
            }
            case (11):
            {
                return {{3, 7}};
                break;
            }
            default:
                MORIS_ASSERT(0, "Invalid edge ordinal specified");
                return moris::Matrix<moris::IndexMat>(0, 0);
                break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_geometric_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Cell_Info_Hex8::get_vertex_loc_coord(moris_index const &aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
            case 0:
            {
                return {{-1.0, -1.0, -1.0}};
                break;
            }
            case 1:
            {
                return {{+1.0, -1.0, -1.0}};
                break;
            }
            case 2:
            {
                return {{+1.0, +1.0, -1.0}};
                break;
            }
            case 3:
            {
                return {{-1.0, +1.0, -1.0}};
                break;
            }
            case 4:
            {
                return {{-1.0, -1.0, +1.0}};
                break;
            }
            case 5:
            {
                return {{+1.0, -1.0, +1.0}};
                break;
            }
            case 6:
            {
                return {{+1.0, +1.0, +1.0}};
                break;
            }
            case 7:
            {
                return {{-1.0, +1.0, +1.0}};
                break;
            }
            default:
            {
                MORIS_ERROR(0, "Invalid vertex ordinal specified");
                return moris::Matrix<moris::DDRMat>(0, 0);
                break;
            }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex8::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
            case (0):
            {
                return {{0, 1}, {1, 5}};
                break;
            }
            case (1):
            {
                return {{1, 2}, {2, 6}};
                break;
            }
            case (2):
            {
                return {{2, 3}, {3, 7}};
                break;
            }
            case (3):
            {
                return {{7, 3}, {3, 0}};
                break;
            }
            case (4):
            {
                return {{0, 3}, {3, 2}};
                break;
            }
            case (5):
            {
                return {{4, 5}, {5, 6}};
                break;
            }
            default:
                MORIS_ERROR(0, "Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0, 0);
                break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex8::compute_cell_size_special(moris::mtk::Cell const *aCell) const
        {
            moris::Cell<Vertex *> tVertices = aCell->get_vertex_pointers();

            Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
            Matrix<DDRMat> tNode6Coords = tVertices(6)->get_coords();

            // FIXME: only works for rectangular cells
            real tLx = std::abs(tNode0Coords(0) - tNode6Coords(0));
            real tLy = std::abs(tNode0Coords(1) - tNode6Coords(1));
            real tLz = std::abs(tNode0Coords(2) - tNode6Coords(2));

            return tLx * tLy * tLz;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex8::compute_cell_size_straight( moris::mtk::Cell const * aCell ) const
        {
            // FIXME: Not consistent with each vertex. Depending what corners are used to be bisected, this won't be consistent

            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();
            Matrix<DDRMat> tNodeCoords(8,3);

            for (uint i = 0; i<8; ++i)
            {
                tNodeCoords( {i,i}, {0,2} ) = tVertices(i)->get_coords()({0,0},{0,2});
            }

            // permutation matrix of vertex IDs to stipulate individual tet4 calculations
            moris::Matrix<DDUMat> tPermutMap = {{1, 0, 5, 2},
                                                {4, 0, 7, 5},
                                                {3, 0, 2, 7},
                                                {6, 2, 5, 7},
                                                {7, 0, 2, 5}};

            // init volume and tet edge vectors
            moris::real tVolume = 0.0;
            Matrix<DDRMat> tNodeCoords10(1,3);
            Matrix<DDRMat> tNodeCoords20(1,3);
            Matrix<DDRMat> tNodeCoords30(1,3);

            // A hex can be broken into 5 separate tets
            for (uint iTet = 0; iTet<5; ++iTet)
            {
                //Assigning Vectors
                tNodeCoords10 = tNodeCoords( { tPermutMap(iTet,1), tPermutMap(iTet,1) }, {0,2} ) - 
                                tNodeCoords( { tPermutMap(iTet,0), tPermutMap(iTet,0) }, {0,2} );
                tNodeCoords20 = tNodeCoords( { tPermutMap(iTet,2), tPermutMap(iTet,2) }, {0,2} ) - 
                                tNodeCoords( { tPermutMap(iTet,0), tPermutMap(iTet,0) }, {0,2} );
                tNodeCoords30 = tNodeCoords( { tPermutMap(iTet,3), tPermutMap(iTet,3) }, {0,2} ) - 
                                tNodeCoords( { tPermutMap(iTet,0), tPermutMap(iTet,0) }, {0,2} );
                
                MORIS_ASSERT(1.0 / 6.0 * dot(tNodeCoords10, cross(tNodeCoords20, tNodeCoords30)) > 0,
                            "Cell_Info_Hex8::compute_cell_size_straight - Determined interior tet "
                            "volume is <=0, suggesting poorly defined nodal coordinates.");

                tVolume +=  1.0 / 6.0 * dot(tNodeCoords10, cross(tNodeCoords20, tNodeCoords30));

            }

            return tVolume;
        }

        // ----------------------------------------------------------------------------------
        moris::real
        Cell_Info_Hex8::compute_cell_side_size(moris::mtk::Cell const *aCell,
                                               moris_index const &aSideOrd) const
        {
            moris::Cell<mtk::Vertex const *> tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            // FIXME: only works for rectangular cells
            Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords1 = tVertices(1)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(3)->get_coords();

            return norm(cross(tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0));
        }

        // ----------------------------------------------------------------------------------

        moris::uint
        Cell_Info_Hex8::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
            case (0):
            {
                return 2;
                break;
            }
            case (1):
            {
                return 3;
                break;
            }
            case (2):
            {
                return 0;
                break;
            }
            case (3):
            {
                return 1;
                break;
            }
            case (4):
            {
                return 5;
                break;
            }
            case (5):
            {
                return 4;
                break;
            }
            default:
            {
                MORIS_ERROR(0, "Invalid side ordinal specified");
                return MORIS_UINT_MAX;
                break;
            }
            }
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Hex8::eval_N(
            const Matrix<DDRMat> &aXi,
            Matrix<DDRMat> &aNXi) const
        {
            // make sure that input is correct
            MORIS_ASSERT(aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size.");

            // unpack xi and eta from input vector
            moris::real xi = aXi(0);
            moris::real eta = aXi(1);
            moris::real zeta = aXi(2);

            // populate output matrix
            aNXi.set_size(1, 8);
            aNXi(0) = -(eta - 1.0) * (xi - 1.0) * (zeta - 1.0) * 0.125;
            aNXi(1) = (eta - 1.0) * (xi + 1.0) * (zeta - 1.0) * 0.125;
            aNXi(2) = -(eta + 1.0) * (xi + 1.0) * (zeta - 1.0) * 0.125;
            aNXi(3) = (eta + 1.0) * (xi - 1.0) * (zeta - 1.0) * 0.125;
            aNXi(4) = (eta - 1.0) * (xi - 1.0) * (zeta + 1.0) * 0.125;
            aNXi(5) = -(eta - 1.0) * (xi + 1.0) * (zeta + 1.0) * 0.125;
            aNXi(6) = (eta + 1.0) * (xi + 1.0) * (zeta + 1.0) * 0.125;
            aNXi(7) = -(eta + 1.0) * (xi - 1.0) * (zeta + 1.0) * 0.125;
        }
    } // namespace mtk
} // namespace moris