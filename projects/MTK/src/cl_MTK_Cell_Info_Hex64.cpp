/*
 * cl_MTK_Cell_Info_Hex64.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */


#include "cl_MTK_Cell_Info_Hex64.hpp"
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
        Cell_Info_Hex64::get_cell_geometry() const
        {
            return Geometry_Type::HEX;
        }

        // ----------------------------------------------------------------------------------
        
        enum CellTopology
        Cell_Info_Hex64::get_cell_topology() const
        {
            return CellTopology::HEX64;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Hex64::get_cell_interpolation_order() const
        {
            return Interpolation_Order::CUBIC;
        }

        // ----------------------------------------------------------------------------------

        enum Integration_Order
        Cell_Info_Hex64::get_cell_integration_order() const
        {
            return Integration_Order::HEX_5x5x5;
        }

        // ----------------------------------------------------------------------------------

        enum CellShape
        Cell_Info_Hex64::compute_cell_shape(moris::mtk::Cell const *aCell) const
        {
            // getting vertices and storing them in a local matrix, since each node will be used a few times
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            // init cell shape
            CellShape tCellShape = CellShape::RECTANGULAR;

            // error threshold
            real tEpsilon = 1.0E-8;

            // init cell of face normals
            moris::Cell< moris::Matrix< DDRMat > > tFaceNormals( 6 );

            // looping through each face
            for ( uint iFace = 0; iFace < 6; iFace++)
            {
                // getting nodes on the face
                moris::Matrix<moris::IndexMat> tFaceNodes = this->get_node_to_face_map( iFace );

                // get the vertex coords
                moris::Matrix< DDRMat > tVertex00 = tVertices( tFaceNodes( 0  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex01 = tVertices( tFaceNodes( 1  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex02 = tVertices( tFaceNodes( 2  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex03 = tVertices( tFaceNodes( 3  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex04 = tVertices( tFaceNodes( 4  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex05 = tVertices( tFaceNodes( 5  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex06 = tVertices( tFaceNodes( 6  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex07 = tVertices( tFaceNodes( 7  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex08 = tVertices( tFaceNodes( 8  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex09 = tVertices( tFaceNodes( 9  ) )->get_coords();
                moris::Matrix< DDRMat > tVertex10 = tVertices( tFaceNodes( 10 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex11 = tVertices( tFaceNodes( 11 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex12 = tVertices( tFaceNodes( 12 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex13 = tVertices( tFaceNodes( 13 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex14 = tVertices( tFaceNodes( 14 ) )->get_coords();
                moris::Matrix< DDRMat > tVertex15 = tVertices( tFaceNodes( 15 ) )->get_coords();

                // get edges to define check plane
                auto tEdge00 = tVertex01 - tVertex00;
                auto tEdge01 = tVertex02 - tVertex01;
                tFaceNormals( iFace ) = cross( tEdge00,tEdge01 );

                // these are the other plane normals produced by the other 6 nodes.
                // keeping node 0 to ensure the plane is at the same offset.
                auto tEdge02 = tVertex03 - tVertex00;
                auto tEdge03 = tVertex04 - tVertex03;
                auto tFaceVec0203 = cross(tEdge02,tEdge03);

                auto tEdge04 = tVertex05 - tVertex00;
                auto tEdge05 = tVertex06 - tVertex05;
                auto tFaceVec0405 = cross(tEdge04,tEdge05);

                auto tEdge06 = tVertex07 - tVertex00;
                auto tEdge07 = tVertex08 - tVertex07;
                auto tFaceVec0607 = cross(tEdge06,tEdge07);

                auto tEdge08 = tVertex09 - tVertex00;
                auto tEdge09 = tVertex10 - tVertex09;
                auto tFaceVec0809 = cross(tEdge08,tEdge09);

                auto tEdge10 = tVertex11 - tVertex00;
                auto tEdge11 = tVertex12 - tVertex11;
                auto tFaceVec1011 = cross(tEdge10,tEdge11);

                auto tEdge12 = tVertex13 - tVertex00;
                auto tEdge13 = tVertex14 - tVertex13;
                auto tFaceVec1213 = cross(tEdge12,tEdge13);

                auto tEdge14 = tVertex15 - tVertex00;
                auto tEdge15 = tVertex01 - tVertex15;
                auto tFaceVec1415 = cross(tEdge14,tEdge15);

                // All three of the plane normals must be parallel in order to be considered straight shape
                if( norm( cross( tFaceNormals( iFace ), tFaceVec0203) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec0405) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec0607) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec0809) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec1011) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec1213) ) > tEpsilon ||
                    norm( cross( tFaceNormals( iFace ), tFaceVec1415) ) > tEpsilon )
                {
                    tCellShape = CellShape::GENERAL;
                    break;
                }
            }

            // since the faces are planar, checking if it is rectangular (aligned), parallel, or straight

            // if the cell shape is straight, determine if it is a parallel cell
            if( tCellShape == CellShape::RECTANGULAR )
            {
                // create the cross products of opposite faces
                auto tCross02 = cross( tFaceNormals(0), tFaceNormals(2) );
                auto tCross13 = cross( tFaceNormals(1), tFaceNormals(3) );
                auto tCross45 = cross( tFaceNormals(4), tFaceNormals(5) );

                // checking if these face normals are parallel
                if( norm( tCross02 ) > tEpsilon ||
                    norm( tCross13 ) > tEpsilon ||
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
                    if ( std::abs(dot( tFaceNormals(0), tFaceNormals(1) )) > tEpsilon ||
                            std::abs(dot( tFaceNormals(0), tFaceNormals(4) )) > tEpsilon ||
                            std::abs(dot( tFaceNormals(1), tFaceNormals(4) )) > tEpsilon ||
                            std::abs( tEdge0(1) ) > tEpsilon ||
                            std::abs( tEdge0(2) ) > tEpsilon ||
                            std::abs( tEdge1(0) ) > tEpsilon ||
                            std::abs( tEdge1(2) ) > tEpsilon ||
                            tEdge0(0)   < 0.0 ||
                            tEdge1(1)   < 0.0 )
                    {
                        tCellShape = CellShape::PARALLEL;
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex64::get_num_verts() const
        {
            return 64;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex64::get_num_facets() const
        {
            return 6;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex64::get_num_edges() const
        {
            return 12;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex64::get_num_verts_per_facet() const
        {
            return 16;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Hex64::get_loc_coord_dim() const
        {
            return 3;
        }

        //-----------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_to_face_map() const
        {
            return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 },
                { 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 },
                { 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 },
                { 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 },
                { 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 },
                { 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 }}; break; }
                case(1):{ return {{ 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 }}; break; }
                case(2):{ return {{ 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 }}; break; }
                case(3):{ return {{ 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 }}; break; }
                case(4):{ return {{ 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 }}; break; }
                case(5):{ return {{ 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }}; break; }
                default:
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
            }
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_to_edge_map() const
        {
            return {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
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
        Cell_Info_Hex64::get_node_to_facet_map() const
        {
            return this->get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_geometric_node_to_facet_map() const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            Cell_Info_Hex8 tHex8;
            return tHex8.get_node_to_face_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Hex64::get_node_map_outward_normal(moris::uint aSideOrdinal)  const
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

        moris::uint
        Cell_Info_Hex64::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return 2; break; }
                case(1):{ return 3; break; }
                case(2):{ return 0; break; }
                case(3):{ return 1; break; }
                case(4):{ return 5; break; }
                case(5):{ return 4; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return MORIS_UINT_MAX;
                    break;
                }
            }
        }

        //-----------------------------------------------------------------------------

        Matrix<DDRMat>
        Cell_Info_Hex64::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
        {
            switch ( aVertexOrdinal )
            {
                case  0: { return {{-1.000000000000000e+00,  -1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case  1: { return {{+1.000000000000000e+00,  -1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case  2: { return {{+1.000000000000000e+00,  +1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case  3: { return {{-1.000000000000000e+00,  +1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case  4: { return {{-1.000000000000000e+00,  -1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case  5: { return {{+1.000000000000000e+00,  -1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case  6: { return {{+1.000000000000000e+00,  +1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case  7: { return {{-1.000000000000000e+00,  +1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case  8: { return {{-3.333333333333334e-01,  -1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case  9: { return {{+3.333333333333333e-01,  -1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case 10: { return {{-1.000000000000000e+00,  -3.333333333333334e-01,  -1.000000000000000e+00}}; break;}
                case 11: { return {{-1.000000000000000e+00,  +3.333333333333333e-01,  -1.000000000000000e+00}}; break;}
                case 12: { return {{-1.000000000000000e+00,  -1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 13: { return {{-1.000000000000000e+00,  -1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 14: { return {{+1.000000000000000e+00,  -3.333333333333334e-01,  -1.000000000000000e+00}}; break;}
                case 15: { return {{+1.000000000000000e+00,  +3.333333333333333e-01,  -1.000000000000000e+00}}; break;}
                case 16: { return {{+1.000000000000000e+00,  -1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 17: { return {{+1.000000000000000e+00,  -1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 18: { return {{+3.333333333333333e-01,  +1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case 19: { return {{-3.333333333333334e-01,  +1.000000000000000e+00,  -1.000000000000000e+00}}; break;}
                case 20: { return {{+1.000000000000000e+00,  +1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 21: { return {{+1.000000000000000e+00,  +1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 22: { return {{-1.000000000000000e+00,  +1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 23: { return {{-1.000000000000000e+00,  +1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 24: { return {{-3.333333333333334e-01,  -1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case 25: { return {{+3.333333333333333e-01,  -1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case 26: { return {{-1.000000000000000e+00,  -3.333333333333334e-01,  +1.000000000000000e+00}}; break;}
                case 27: { return {{-1.000000000000000e+00,  +3.333333333333333e-01,  +1.000000000000000e+00}}; break;}
                case 28: { return {{+1.000000000000000e+00,  -3.333333333333334e-01,  +1.000000000000000e+00}}; break;}
                case 29: { return {{+1.000000000000000e+00,  +3.333333333333333e-01,  +1.000000000000000e+00}}; break;}
                case 30: { return {{+3.333333333333333e-01,  +1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case 31: { return {{-3.333333333333334e-01,  +1.000000000000000e+00,  +1.000000000000000e+00}}; break;}
                case 32: { return {{-3.333333333333334e-01,  -3.333333333333334e-01,  -1.000000000000000e+00}}; break;}
                case 33: { return {{-3.333333333333334e-01,  +3.333333333333333e-01,  -1.000000000000000e+00}}; break;}
                case 34: { return {{+3.333333333333333e-01,  +3.333333333333333e-01,  -1.000000000000000e+00}}; break;}
                case 35: { return {{+3.333333333333333e-01,  -3.333333333333334e-01,  -1.000000000000000e+00}}; break;}
                case 36: { return {{-3.333333333333334e-01,  -1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 37: { return {{+3.333333333333333e-01,  -1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 38: { return {{+3.333333333333333e-01,  -1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 39: { return {{-3.333333333333334e-01,  -1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 40: { return {{-1.000000000000000e+00,  -3.333333333333334e-01,  -3.333333333333334e-01}}; break;}
                case 41: { return {{-1.000000000000000e+00,  -3.333333333333334e-01,  +3.333333333333333e-01}}; break;}
                case 42: { return {{-1.000000000000000e+00,  +3.333333333333333e-01,  +3.333333333333333e-01}}; break;}
                case 43: { return {{-1.000000000000000e+00,  +3.333333333333333e-01,  -3.333333333333334e-01}}; break;}
                case 44: { return {{+1.000000000000000e+00,  -3.333333333333334e-01,  -3.333333333333334e-01}}; break;}
                case 45: { return {{+1.000000000000000e+00,  +3.333333333333333e-01,  -3.333333333333334e-01}}; break;}
                case 46: { return {{+1.000000000000000e+00,  +3.333333333333333e-01,  +3.333333333333333e-01}}; break;}
                case 47: { return {{+1.000000000000000e+00,  -3.333333333333334e-01,  +3.333333333333333e-01}}; break;}
                case 48: { return {{+3.333333333333333e-01,  +1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 49: { return {{-3.333333333333334e-01,  +1.000000000000000e+00,  -3.333333333333334e-01}}; break;}
                case 50: { return {{-3.333333333333334e-01,  +1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 51: { return {{+3.333333333333333e-01,  +1.000000000000000e+00,  +3.333333333333333e-01}}; break;}
                case 52: { return {{-3.333333333333334e-01,  -3.333333333333334e-01,  +1.000000000000000e+00}}; break;}
                case 53: { return {{+3.333333333333333e-01,  -3.333333333333334e-01,  +1.000000000000000e+00}}; break;}
                case 54: { return {{+3.333333333333333e-01,  +3.333333333333333e-01,  +1.000000000000000e+00}}; break;}
                case 55: { return {{-3.333333333333334e-01,  +3.333333333333333e-01,  +1.000000000000000e+00}}; break;}
                case 56: { return {{-3.333333333333334e-01,  -3.333333333333334e-01,  -3.333333333333334e-01}}; break;}
                case 57: { return {{+3.333333333333333e-01,  -3.333333333333334e-01,  -3.333333333333334e-01}}; break;}
                case 58: { return {{+3.333333333333333e-01,  +3.333333333333333e-01,  -3.333333333333334e-01}}; break;}
                case 59: { return {{-3.333333333333334e-01,  +3.333333333333333e-01,  -3.333333333333334e-01}}; break;}
                case 60: { return {{-3.333333333333334e-01,  -3.333333333333334e-01,  +3.333333333333333e-01}}; break;}
                case 61: { return {{+3.333333333333333e-01,  -3.333333333333334e-01,  +3.333333333333333e-01}}; break;}
                case 62: { return {{+3.333333333333333e-01,  +3.333333333333333e-01,  +3.333333333333333e-01}}; break;}
                case 63: { return {{-3.333333333333334e-01,  +3.333333333333333e-01,  +3.333333333333333e-01}}; break;}
                default:
                {
                    MORIS_ERROR(0,"Invalid vertex ordinal specified");
                    return moris::Matrix<moris::DDRMat>(0,0);
                    break;
                }
            }
        }

        //-----------------------------------------------------------------------------
        moris::real
        Cell_Info_Hex64::compute_cell_size_special( moris::mtk::Cell const * aCell ) const
        {
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
            Matrix<DDRMat> tNode6Coords = tVertices(6)->get_coords();

            // FIXME: only works for rectangular cells
            real tLx = std::abs(tNode0Coords(0) - tNode6Coords(0));
            real tLy = std::abs(tNode0Coords(1) - tNode6Coords(1));
            real tLz = std::abs(tNode0Coords(2) - tNode6Coords(2));

            return tLx*tLy*tLz;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Hex64::compute_cell_side_size(
                moris::mtk::Cell const * aCell ,
                moris_index      const & aSideOrd) const
        {
            moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            // FIXME: only works for rectangular cells
            Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
            Matrix<DDRMat> tNodeCoords1 = tVertices(1)->get_coords();
            Matrix<DDRMat> tNodeCoords2 = tVertices(3)->get_coords();

            return norm( cross( tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0 ) );
        }

        // ----------------------------------------------------------------------------------

        void
        Cell_Info_Hex64::eval_N( const Matrix< DDRMat > & aXi,
                Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 3, "HEX64 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta and zeta from input vector
            real   xi = aXi( 0 );
            real  eta = aXi( 1 );
            real zeta = aXi( 2 );

            real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

            real c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
            real c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
            real c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
            real c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

            // populate matrix with values
            aNXi.set_size(1,64);
            aNXi(  0 ) = a0 * b0 * c0;
            aNXi(  1 ) = a3 * b0 * c0;
            aNXi(  2 ) = a3 * b3 * c0;
            aNXi(  3 ) = a0 * b3 * c0;
            aNXi(  4 ) = a0 * b0 * c3;
            aNXi(  5 ) = a3 * b0 * c3;
            aNXi(  6 ) = a3 * b3 * c3;
            aNXi(  7 ) = a0 * b3 * c3;
            aNXi(  8 ) = a1 * b0 * c0;
            aNXi(  9 ) = a2 * b0 * c0;
            aNXi( 10 ) = a0 * b1 * c0;
            aNXi( 11 ) = a0 * b2 * c0;
            aNXi( 12 ) = a0 * b0 * c1;
            aNXi( 13 ) = a0 * b0 * c2;
            aNXi( 14 ) = a3 * b1 * c0;
            aNXi( 15 ) = a3 * b2 * c0;
            aNXi( 16 ) = a3 * b0 * c1;
            aNXi( 17 ) = a3 * b0 * c2;
            aNXi( 18 ) = a2 * b3 * c0;
            aNXi( 19 ) = a1 * b3 * c0;
            aNXi( 20 ) = a3 * b3 * c1;
            aNXi( 21 ) = a3 * b3 * c2;
            aNXi( 22 ) = a0 * b3 * c1;
            aNXi( 23 ) = a0 * b3 * c2;
            aNXi( 24 ) = a1 * b0 * c3;
            aNXi( 25 ) = a2 * b0 * c3;
            aNXi( 26 ) = a0 * b1 * c3;
            aNXi( 27 ) = a0 * b2 * c3;
            aNXi( 28 ) = a3 * b1 * c3;
            aNXi( 29 ) = a3 * b2 * c3;
            aNXi( 30 ) = a2 * b3 * c3;
            aNXi( 31 ) = a1 * b3 * c3;
            aNXi( 32 ) = a1 * b1 * c0;
            aNXi( 33 ) = a1 * b2 * c0;
            aNXi( 34 ) = a2 * b2 * c0;
            aNXi( 35 ) = a2 * b1 * c0;
            aNXi( 36 ) = a1 * b0 * c1;
            aNXi( 37 ) = a2 * b0 * c1;
            aNXi( 38 ) = a2 * b0 * c2;
            aNXi( 39 ) = a1 * b0 * c2;
            aNXi( 40 ) = a0 * b1 * c1;
            aNXi( 41 ) = a0 * b1 * c2;
            aNXi( 42 ) = a0 * b2 * c2;
            aNXi( 43 ) = a0 * b2 * c1;
            aNXi( 44 ) = a3 * b1 * c1;
            aNXi( 45 ) = a3 * b2 * c1;
            aNXi( 46 ) = a3 * b2 * c2;
            aNXi( 47 ) = a3 * b1 * c2;
            aNXi( 48 ) = a2 * b3 * c1;
            aNXi( 49 ) = a1 * b3 * c1;
            aNXi( 50 ) = a1 * b3 * c2;
            aNXi( 51 ) = a2 * b3 * c2;
            aNXi( 52 ) = a1 * b1 * c3;
            aNXi( 53 ) = a2 * b1 * c3;
            aNXi( 54 ) = a2 * b2 * c3;
            aNXi( 55 ) = a1 * b2 * c3;
            aNXi( 56 ) = a1 * b1 * c1;
            aNXi( 57 ) = a2 * b1 * c1;
            aNXi( 58 ) = a2 * b2 * c1;
            aNXi( 59 ) = a1 * b2 * c1;
            aNXi( 60 ) = a1 * b1 * c2;
            aNXi( 61 ) = a2 * b1 * c2;
            aNXi( 62 ) = a2 * b2 * c2;
            aNXi( 63 ) = a1 * b2 * c2;
        }
    }
}
