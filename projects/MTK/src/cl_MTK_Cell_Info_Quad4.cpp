/*
 * cl_MTK_Cell_Info_Quad4.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------

        enum Geometry_Type
        Cell_Info_Quad4::get_cell_geometry() const
        {
            return Geometry_Type::QUAD;
        }

        // ----------------------------------------------------------------------------------

        enum CellTopology
        Cell_Info_Quad4::get_cell_topology() const
        {
            return CellTopology::QUAD4;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Quad4::get_cell_interpolation_order() const
        {
            return Interpolation_Order::LINEAR;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Quad4::get_cell_integration_order() const
        {
            return Integration_Order::QUAD_2x2;
        }

        // ----------------------------------------------------------------------------------

        enum CellTopology
        Cell_Info_Quad4::get_cell_topology() const
        {
            return CellTopology::QUAD4;
        }

        //-----------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Quad4::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // error threshold
            real tEpsilon = 1.0E-8;

            // get the vertex coordinates
            moris::Matrix< DDRMat > tVertex0 = tVertices( 0 )->get_coords();
            moris::Matrix< DDRMat > tVertex1 = tVertices( 1 )->get_coords();
            moris::Matrix< DDRMat > tVertex2 = tVertices( 2 )->get_coords();
            moris::Matrix< DDRMat > tVertex3 = tVertices( 3 )->get_coords();

            // getting edge vectors
            moris::Matrix< DDRMat > tEdge0 = tVertex1 - tVertex0;
            moris::Matrix< DDRMat > tEdge1 = tVertex2 - tVertex1;
            moris::Matrix< DDRMat > tEdge2 = tVertex3 - tVertex2;
            moris::Matrix< DDRMat > tEdge3 = tVertex0 - tVertex3;

            // cross products of opposite edges
            auto tCross02 = tEdge0(0)*tEdge2(1)-tEdge0(1)*tEdge2(0);
            auto tCross13 = tEdge1(0)*tEdge3(1)-tEdge1(1)*tEdge3(0);

            // check if opposite edges are parallel
            if ( std::abs( tCross02 ) > tEpsilon ||
                 std::abs( tCross13 ) > tEpsilon )
            {
                return CellShape::STRAIGHT;
            }

            // parallelogram cell
            else
            {
                // if edge 1 is parallel to the x axis and perpindicular to the adjacent edge
                if( std::abs( tEdge0(1) )          < tEpsilon &&
                    std::abs(dot( tEdge0,tEdge1 )) < tEpsilon )
                {
                    return CellShape::RECTANGULAR;
                }
                else
                {
                    return CellShape::PARALLEL;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad4::get_num_verts() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad4::get_num_facets() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad4::get_num_edges() const
        {
            return 4;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad4::get_num_verts_per_facet() const
        {
            return 2;
        }
        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Quad4::get_loc_coord_dim() const
        {
            return 2;
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_to_face_map() const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_to_edge_map() const
        {
            return {{0,1}, {1,2}, {2,3}, {3,0}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_to_facet_map() const
        {
            return this->get_node_to_edge_map();
        }

        //-----------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_geometric_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_facet_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------
        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
                case(0):{ return {{0, 1}}; break; }
                case(1):{ return {{1, 2}}; break; }
                case(2):{ return {{2, 3}}; break; }
                case(3):{ return {{3, 0}}; break; }
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
        Cell_Info_Quad4::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Quad4::get_node_map_outward_normal(moris::uint aSideOrdinal) const
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
        Cell_Info_Quad4::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
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

        moris::Cell< moris_index >
        Cell_Info_Quad4::get_vertex_path_to_entity_rank_and_ordinal(
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
                                
                                case 1: return { 0, 1 };
                                case 2: return { 0, 3 };
                                case 3: return { 3, 1 };

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }

                        // node to edge paths
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case  0: return { 0, 1 };
                                case  1: return { 0, 3 };
                                case  2: return { 0, 3 };
                                case  3: return { 3, 1 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }
                        
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch aOtherEntityRank
                } // end: case aVertexOrdinal == 0
                case 1:
                {
                    switch ( aOtherEntityRank )
                    {
                        // node to node paths
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 1:
                                {
                                    MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                    return 0;
                                }
                                
                                case 0: return { 0, 1 };
                                case 2: return { 1, 1 };
                                case 3: return { 0, 3 }; 

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }

                        // node to edge paths
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case  0: return { 0, 1 };
                                case  1: return { 1, 1 };
                                case  2: return { 0, 3 };
                                case  3: return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }
                        
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch aOtherEntityRank
                } // end: case aVertexOrdinal == 0
                case 2:
                {
                    switch ( aOtherEntityRank )
                    {
                        // node to node paths
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 2:
                                {
                                    MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                    return 0;
                                }
                                
                                case 0: return { 0, 3 };
                                case 1: return { 1, 1 };
                                case 3: return { 2, 1 };

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }

                        // node to edge paths
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case  0: return { 0, 3 };
                                case  1: return { 1, 1 };
                                case  2: return { 2, 1 };
                                case  3: return { 0, 3 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }
                        
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch aOtherEntityRank
                } // end: case aVertexOrdinal == 0
                case 3:
                {
                    switch ( aOtherEntityRank )
                    {
                        // node to node paths
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 3:
                                {
                                    MORIS_ERROR( 0, "No Path between a vertex and itself" );
                                    return 0;
                                }
                                
                                case 0: return { 3, 1 };
                                case 1: return { 0, 3 };
                                case 2: return { 2, 1 };

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other node ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }

                        // node to edge paths
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case  0: return { 0, 3 };
                                case  1: return { 0, 3 };
                                case  2: return { 2, 1 };
                                case  3: return { 3, 1 };
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for hex8" );
                                    return moris::Cell< moris_index >( 0 );
                                }
                            }
                        }
                        
                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for hex8" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch aOtherEntityRank
                } // end: case aVertexOrdinal == 0

                default:
                {
                    MORIS_ERROR( 0, "Invalid vertex ordinal for hex8" );
                    return moris::Cell< moris_index >( 0 );
                }
            } // end: switch aVertexOrdinal
        }

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index >
        Cell_Info_Quad4::get_edge_path_to_entity_rank_and_ordinal(
            moris_index aEdgeOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
        {
            switch ( aEdgeOrdinal )
            {
                // edge 0
                case 0:
                {
                    switch ( aOtherEntityRank )
                    {
                        // other entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return { 0, 1 };
                                case 1: return { 0, 3 };
                                case 2: return { 0, 3 };
                                case 3: return { 0, 3 };
                                
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return moris::Cell< moris_index >( 0 );
                                }   
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 - edge has only path to edges" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 0
                
                // edge 1
                case 1:
                {
                    switch ( aOtherEntityRank )
                    {
                        // other entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return { 0, 3 };
                                case 1: return { 0, 1 };
                                case 2: return { 0, 3 };
                                case 3: return { 0, 3 };
                                
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return moris::Cell< moris_index >( 0 );
                                }   
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 - edge only has path to edges" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 1

                // edge 2
                case 2:
                {
                    switch ( aOtherEntityRank )
                    {
                        // other entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return { 0, 3 };
                                case 1: return { 0, 3 };
                                case 2: return { 0, 1 };
                                case 3: return { 0, 3 };
                                
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return moris::Cell< moris_index >( 0 );
                                }   
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 - edge only has path to edges" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 2

                // edge 3
                case 3:
                {
                    switch ( aOtherEntityRank )
                    {
                        // other entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return { 0, 3 };
                                case 1: return { 0, 3 };
                                case 2: return { 0, 3 };
                                case 3: return { 0, 1 };
                                
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return moris::Cell< moris_index >( 0 );
                                }   
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 - edge only has path to edges" );
                            return moris::Cell< moris_index >( 0 );
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 3

                default:
                {
                    MORIS_ERROR( 0, "Invalid vertex ordinal for QUAD4" );
                    return moris::Cell< moris_index >( 0 );
                }
            }
        }

        // ----------------------------------------------------------------------------------

        bool
        Cell_Info_Quad4::is_entity_connected_to_facet(
            moris_index aFacetOrdinal,
            moris_index aOtherEntityOrdinal,
            moris_index aOtherEntityRank ) const
        {
            switch ( aFacetOrdinal )
            {
                // Edge Ordinal 0
                case 0:
                {
                    switch ( aOtherEntityRank )
                    {
                        // Other Entity is Vertex
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return true;
                                case 1: return true;
                                case 2: return false;
                                case 3: return false;

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other vertex ordinal for Quad4" );
                                    return false;
                                }
                            }
                        }

                        // Other Entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return true;
                                case 1: return false;
                                case 2: return false;
                                case 3: return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return false;
                                }
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 (must be 0-Vertex or 1-Edge" );
                            return false;
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 0

                // Edge Ordinal 1
                case 1:
                {
                    switch ( aOtherEntityRank )
                    {
                        // Other Entity is Vertex
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return false;
                                case 1: return true;
                                case 2: return true;
                                case 3: return false;

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other vertex ordinal for Quad4" );
                                    return false;
                                }
                            }
                        }

                        // Other Entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return false;
                                case 1: return true;
                                case 2: return false;
                                case 3: return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return false;
                                }
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 (must be 0-Vertex or 1-Edge" );
                            return false;
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 1

                // Edge Ordinal 2
                case 2:
                {
                    switch ( aOtherEntityRank )
                    {
                        // Other Entity is Vertex
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return false;
                                case 1: return false;
                                case 2: return true;
                                case 3: return true;

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other vertex ordinal for Quad4" );
                                    return false;
                                }
                            }
                        }

                        // Other Entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return false;
                                case 1: return false;
                                case 2: return true;
                                case 3: return false;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return false;
                                }
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 (must be 0-Vertex or 1-Edge" );
                            return false;
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 2

                // Edge Ordinal 3
                case 3:
                {
                    switch ( aOtherEntityRank )
                    {
                        // Other Entity is Vertex
                        case 0:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return true;
                                case 1: return false;
                                case 2: return false;
                                case 3: return true;

                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other vertex ordinal for Quad4" );
                                    return false;
                                }
                            }
                        }

                        // Other Entity is Edge
                        case 1:
                        {
                            switch ( aOtherEntityOrdinal )
                            {
                                case 0: return false;
                                case 1: return false;
                                case 2: return false;
                                case 3: return true;
                                default:
                                {
                                    MORIS_ERROR( 0, "Invalid other edge ordinal for QUAD4" );
                                    return false;
                                }
                            }
                        }

                        default:
                        {
                            MORIS_ERROR( 0, "Invalid other entity rank for QUAD4 (must be 0-Vertex or 1-Edge" );
                            return false;
                        }
                    } // end: switch: aOtherEntityRank
                } // end: case: edge ordinal 3

                default:
                {
                    MORIS_ERROR( 0, "Invalid facet ordinal for QUAD4" );
                    return false;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Cell_Info_Quad4::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch (aVertexOrdinal)
            {
                case 0: { return {{-1.000000000000000e+00,  -1.000000000000000e+00}}; break; }
                case 1: { return {{+1.000000000000000e+00,  -1.000000000000000e+00}}; break; }
                case 2: { return {{+1.000000000000000e+00,  +1.000000000000000e+00}}; break; }
                case 3: { return {{-1.000000000000000e+00,  +1.000000000000000e+00}}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid vertex ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad4::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            const Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            const Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2,"Cell_Info_Quad4::compute_cell_size_special only works in 2D.\n");

            // FIXME: only works for rectangular cells
            real tLx = std::abs(tNodeCoords0(0) - tNodeCoords2(0));
            real tLy = std::abs(tNodeCoords0(1) - tNodeCoords2(1));

            real tArea = tLx*tLy;

            return tArea;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad4::compute_cell_size_straight( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            const Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            const Matrix<DDRMat> tNodeCoords1 = tVertices(1)->get_coords();
            const Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();
            const Matrix<DDRMat> tNodeCoords3 = tVertices(3)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2,"Cell_Info_Quad4::compute_cell_size_straight only works in 2D.\n");

            // computes the cross product of the 2 triangles and adds them. some simplifications made.
            real tArea = 0.5 * ( tNodeCoords0(0) *   tNodeCoords1(1) +
                                 tNodeCoords1(0) * ( tNodeCoords2(1) - tNodeCoords0(0) ) - 
                                 tNodeCoords2(0) *   tNodeCoords1(1) +
                                 tNodeCoords2(0) *   tNodeCoords3(1) +
                                 tNodeCoords3(0) * ( tNodeCoords0(1) - tNodeCoords2(1) ) -
                                 tNodeCoords0(0) *   tNodeCoords3(1) );

            return tArea;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad4::compute_cell_size_deriv( moris::mtk::Cell const * aCell, uint aLocalVertexID, uint aDirection ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // permutation vector used to index correct vertices
            moris::Matrix< DDUMat > tVertIndexMap = {{1,2,3,0,1,2}};
            moris::Matrix< DDUMat > tDirIndexMap = {{1,0}};

            // Getting adjacent vertices to vertex of interest
            const Matrix<DDRMat> tNodeCoordsA = tVertices( tVertIndexMap( aLocalVertexID ))->get_coords();
            const Matrix<DDRMat> tNodeCoordsB = tVertices( tVertIndexMap( aLocalVertexID + 2 ))->get_coords();

            MORIS_ASSERT(tNodeCoordsA.numel() == 2,"Cell_Info_Quad4::compute_cell_size_deriv only works in 2D.\n");
            MORIS_ASSERT( aDirection < 2,"Cell_Info_Quad4::compute_cell_size_deriv directions can only be 0 or 1.\n");
            MORIS_ASSERT( aLocalVertexID < 4,"Cell_Info_Quad4::compute_cell_size_deriv vertex IDs must be 0, 1, 2, or 3.\n");

            // computes the derivative of the area wrt to the single dof/direction.
            moris::real tAreaDeriv = 0.5 * std::pow(-1.0, aDirection) *
                                           ( tNodeCoordsA( tDirIndexMap( aDirection ) ) -
                                             tNodeCoordsB( tDirIndexMap( aDirection ) ) );

            return tAreaDeriv;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad4::compute_cell_side_size(
                moris::mtk::Cell const * aCell ,
                moris_index      const & aSideOrd) const
        {
            const moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            return moris::norm(tVertices(1)->get_coords() - tVertices(0)->get_coords());
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Quad4::compute_cell_side_size_deriv(
                moris::mtk::Cell const * aCell,
                moris_index const      & aSideOrd,
                uint                     aLocalVertexID,
                uint                     aDirection ) const
        {
            MORIS_ERROR(false,"compute_cell_side_size_deriv not implemented for Cell_Info_Quad4 yet.");
            return 0.0;
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Quad4::eval_N(
                const Matrix< DDRMat > & aXi,
                Matrix< DDRMat >       & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // populate matrix with values
            aNXi.set_size( 1, 4 );
            aNXi( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            aNXi( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            aNXi( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            aNXi( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
        }
    }
}
