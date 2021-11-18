/*
 * cl_MTK_Cell_Info_Tri3.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
#include "op_times.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------

        enum Geometry_Type
        Cell_Info_Tri3::get_cell_geometry() const
        {
            return Geometry_Type::TRI;
        }

        // ----------------------------------------------------------------------------------

        enum CellTopology
        Cell_Info_Tri3::get_cell_topology() const
        {
            return CellTopology::TRI3;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Tri3::get_cell_interpolation_order() const
        {
            return Interpolation_Order::LINEAR;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Tri3::get_cell_integration_order() const
        {
            return Integration_Order::TRI_7;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Tri3::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            return CellShape::STRAIGHT;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_verts() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_facets() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_edges() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_verts_per_facet() const
        {
            return 2;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_loc_coord_dim() const
        {
            return 3;
        }
        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_face_map() const
        {
            MORIS_ERROR(0, "Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0, 0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_edge_map() const
        {
            return {{0, 1}, {1, 2}, {2, 0}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_facet_map() const
        {
            return this->get_node_to_edge_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0, "Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0, 0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
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
                return {{2, 0}};
                break;
            }
            default:
            {
                MORIS_ASSERT(0, "Invalid edge ordinal specified");
                return moris::Matrix<moris::IndexMat>(0, 0);
                break;
            }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
            case (0):
            {
                return {{1, 0}};
                break;
            }
            case (1):
            {
                return {{2, 1}};
                break;
            }
            case (2):
            {
                return {{0, 2}};
                break;
            }
            default:
            {
                MORIS_ERROR(0, "Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0, 0);
                break;
            }
            }
        }

        // ----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Cell_Info_Tri3::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
                case 0: { return {{+1.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00}}; break; }
                case 1: { return {{ 0.000000000000000e+00, +1.000000000000000e+00,  0.000000000000000e+00}}; break; }
                case 2: { return {{ 0.000000000000000e+00,  0.000000000000000e+00, +1.000000000000000e+00}}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid vertex ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real Cell_Info_Tri3::compute_cell_size_special(moris::mtk::Cell const *aCell) const
        {
            // cell coordinates
            moris::Cell<Vertex *> tVertices = aCell->get_vertex_pointers();

            const Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2, "Cell_Info_Tri3::compute_cell_size_special only works in 2D.\n");

            const Matrix<DDRMat> tNodeCoords10 = tVertices(1)->get_coords() - tNodeCoords0;
            const Matrix<DDRMat> tNodeCoords20 = tVertices(2)->get_coords() - tNodeCoords0;

            real tArea = 0.5 *std::abs(tNodeCoords10(0) * tNodeCoords20(1) - tNodeCoords20(0) * tNodeCoords10(1));

            return tArea; 
        }

        // ----------------------------------------------------------------------------------

        moris::real Cell_Info_Tri3::compute_cell_size_straight(moris::mtk::Cell const *aCell) const
        {
            return compute_cell_size_special(aCell);
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Tri3::compute_cell_size_deriv( moris::mtk::Cell const * aCell, uint aLocalVertexID, uint aDirection ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // permutation vector used to index correct vertices
            moris::Matrix< DDUMat > tVertIndexMap = {{1,2,0,1}};
            moris::Matrix< DDUMat > tDirIndexMap = {{1,0}};

            // Getting adjacent vertices to vertex of interest
            const Matrix<DDRMat> tNodeCoordsA = tVertices( tVertIndexMap( aLocalVertexID ))->get_coords();
            const Matrix<DDRMat> tNodeCoordsB = tVertices( tVertIndexMap( aLocalVertexID + 1 ))->get_coords();

            MORIS_ASSERT(tNodeCoordsA.numel() == 2,"Cell_Info_Tri3::compute_cell_size_deriv only works in 2D.\n");
            MORIS_ASSERT( aDirection < 2,"Cell_Info_Tri3::compute_cell_size_deriv directions can only be 0 or 1.\n");
            MORIS_ASSERT( aLocalVertexID < 3,"Cell_Info_Tri3::compute_cell_size_deriv vertex IDs must be 0, 1, or 2.\n");

            // computes the derivative of the area wrt to the single dof/direction.
            moris::real tAreaDeriv = 0.5 * std::pow( -1.0, aDirection ) *
                    ( tNodeCoordsA( tDirIndexMap( aDirection ) ) -
                            tNodeCoordsB( tDirIndexMap( aDirection ) ) );

            return tAreaDeriv;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Tri3::compute_cell_side_size(
            moris::mtk::Cell const *aCell,
            moris_index const &aSideOrd) const
        {
            moris::Cell<mtk::Vertex const *> tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            Matrix<DDRMat> tLVec = tVertices(1)->get_coords() - tVertices(0)->get_coords();

            return moris::norm(tLVec);
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Tri3::compute_cell_side_size_deriv(
                moris::mtk::Cell const * aCell,
                moris_index const      & aSideOrd,
                uint                     aLocalVertexID,
                uint                     aDirection ) const
        {
            moris::Cell<mtk::Vertex const *> tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            // Getting adjacent vertices to vertex of interest
            const Matrix<DDRMat> tNodeCoordsA = tVertices( 0 )->get_coords();
            const Matrix<DDRMat> tNodeCoordsB = tVertices( 1 )->get_coords();

            // Computing the side length
            moris::real tLength = moris::norm(tNodeCoordsB - tNodeCoordsA);

            MORIS_ASSERT(tNodeCoordsA.numel() == 2,"Cell_Info_Tri3::compute_cell_side_size_deriv only works in 2D.\n");
            MORIS_ASSERT( aDirection < 2,"Cell_Info_Tri3::compute_cell_side_size_deriv directions can only be 0 or 1.\n");
            MORIS_ASSERT( aLocalVertexID < 2,"Cell_Info_Tri3::compute_cell_side_size_deriv vertex IDs must be 0, 1.\n");

            // computes the derivative of the area wrt to the single dof/direction.
            moris::real tLengthDeriv = std::pow( -1.0, aLocalVertexID + 1 ) *
                    ( tNodeCoordsB( aDirection ) - tNodeCoordsA( aDirection ) ) / tLength;

            return tLengthDeriv;
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }


    } // namespace mtk
} // namespace moris
