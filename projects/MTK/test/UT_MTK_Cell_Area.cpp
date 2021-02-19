/*
 * UT_MTK_Cell_Area.cpp
 *
 *  Created on: Jan. 20, 2021
 *      Author: gates
 */

#include "catch.hpp"
#include "cl_Communication_Tools.hpp"

// base class
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Cell_STK.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Vertex_STK.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "IP/cl_MTK_Space_Interpolator.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_reindex_mat.hpp"

namespace moris
{
    namespace mtk
    {

        TEST_CASE("MTK Hex8 Rectangular", "[MTK],[Rect_Hex8]")
        {
            if (par_size() <= 1)
            {
                // define a HEX8 space element, i.e. space coordinates xHat
                Matrix<DDRMat> tXHat = {{0.0, 0.0, 0.0},
                                        {1.2, 0.0, 0.0},
                                        {1.2, 1.0, 0.0},
                                        {0.0, 1.0, 0.0},
                                        {0.0, 0.0, 1.5},
                                        {1.2, 0.0, 1.5},
                                        {1.2, 1.0, 1.5},
                                        {0.0, 1.0, 1.5}};

                // the HEX8 interpolation element in space and time param coordinates xiHat
                Matrix<DDRMat> tXiHat = {{-1.0, -1.0, -1.0},
                                         {1.0, -1.0, -1.0},
                                         {1.0, 1.0, -1.0},
                                         {-1.0, 1.0, -1.0},
                                         {-1.0, -1.0, 1.0},
                                         {1.0, -1.0, 1.0},
                                         {1.0, 1.0, 1.0},
                                         {-1.0, 1.0, 1.0}};

                // integration mesh interpolation rule
                Interpolation_Rule tInterpRule(Geometry_Type::HEX,
                                               Interpolation_Type::LAGRANGE,
                                               Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               Interpolation_Order::QUADRATIC);

                // create a space interpolator
                Space_Interpolator tSpaceInterpolator(tInterpRule);

                // set the coefficients xHat, tHat
                tSpaceInterpolator.set_space_coeff(tXHat);

                // set the coefficients xiHat, tauHat
                tSpaceInterpolator.set_space_param_coeff(tXiHat);

                // Interpolation location
                Matrix<DDRMat> tX = {{0.4}, {0.4}, {0.4}};

                // Set the space
                tSpaceInterpolator.set_space(tX);

                // calculating det J and inv J
                real tDetJ = tSpaceInterpolator.space_det_J();
                Matrix<DDRMat> tInvJac = tSpaceInterpolator.inverse_space_jacobian();

                // setting rectangular flag
                bool tRectangular = true;
                tSpaceInterpolator.set_rectangular(tRectangular);
                tSpaceInterpolator.reset_eval_flags();

                // calculating det J and inv J using rectangular eval
                real tDetJ_Rect = tSpaceInterpolator.space_det_J();
                Matrix<DDRMat> tInvJac_Rect = tSpaceInterpolator.inverse_space_jacobian();

                // checking det J calc
                CHECK(tDetJ == Approx(tDetJ_Rect));

                // checking inverse jacobian
                for (uint i = 0; i <= 2; i++)
                {
                    for (uint j = 0; j <= 2; j++)
                    {
                        CHECK(tInvJac(i, j) == Approx(tInvJac_Rect(i, j)));
                    }
                }
            }
        }

        TEST_CASE("MTK Quad4 Rectangular", "[MTK],[Rect_Quad4]")
        {
            if (par_size() <= 1)
            {
                // define a QUAD4 space element, i.e. space coordinates xHat
                Matrix<DDRMat> tXHat = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.5}, {0.0, 1.5}};

                // the QUAD4 interpolation element in space and time param coordinates xiHat
                Matrix<DDRMat> tXiHat = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};

                // integration mesh interpolation rule
                Interpolation_Rule tInterpRule(Geometry_Type::QUAD,
                                               Interpolation_Type::LAGRANGE,
                                               Interpolation_Order::LINEAR,
                                               Interpolation_Type::LAGRANGE,
                                               Interpolation_Order::QUADRATIC);

                // create a space interpolator
                Space_Interpolator tSpaceInterpolator(tInterpRule);

                // set the coefficients xHat, tHat
                tSpaceInterpolator.set_space_coeff(tXHat);

                // set the coefficients xiHat, tauHat
                tSpaceInterpolator.set_space_param_coeff(tXiHat);

                // Interpolation location
                Matrix<DDRMat> tX = {{0.4}, {0.4}};

                // Set the space
                tSpaceInterpolator.set_space(tX);

                // calculating detJ and the inverse jacobian
                real tDetJ = tSpaceInterpolator.space_det_J();
                Matrix<DDRMat> tInvJac = tSpaceInterpolator.inverse_space_jacobian();

                // setting rectangular flag
                bool tRectangular = true;
                tSpaceInterpolator.set_rectangular(tRectangular);
                tSpaceInterpolator.reset_eval_flags();

                // calculating detJ and inv J using rectangular calc
                real tDetJ_Rect = tSpaceInterpolator.space_det_J();
                Matrix<DDRMat> tInvJac_Rect = tSpaceInterpolator.inverse_space_jacobian();

                // checking det j
                CHECK(tDetJ == Approx(tDetJ_Rect));

                // checking inverse jacobian
                for (uint i = 0; i <= 1; i++)
                {
                    for (uint j = 0; j <= 1; j++)
                    {
                        CHECK(tInvJac(i, j) == Approx(tInvJac_Rect(i, j)));
                    }
                }
            }
        }

        TEST_CASE("MTK Cell Area Integration", "[MTK],[Area_Integration],[Area_Quad4]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.0, 0.0},
                                          {1.7, 1.7},
                                          {0.7, 1.5}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                // specify the local to global map
                //Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD4;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                //aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                std::string tFileOutput = "./mtk_quad4_cell_area_ut.exo";
                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // Setup Node Vertices (note: this data structure will be in the STK_Implementation
                moris::Cell<Vertex *> tElementVertices;
                for (size_t i = 0; i < tNodeIndices.numel(); i++)
                {
                    tElementVertices.push_back(new Vertex_STK(tNodeIds(i), tNodeIndices(i), tMesh));
                }

                // Setup cell associated with element index 0
                Cell_Info_Quad4 tQuad4;
                Cell_STK tCell(&tQuad4, 1, 0, tElementVertices, tMesh);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 4);
                Interpolation_Order tInterpOrder = tCell.get_interpolation_order();
                CHECK(tInterpOrder == Interpolation_Order::LINEAR);

                Integration_Order tIntegOrder = tCell.get_integration_order();
                CHECK(tIntegOrder == Integration_Order::QUAD_2x2);

                CHECK(tCell.compute_cell_measure_general() == 1.53);
            }
        }

        TEST_CASE("MTK Cell Area Integration Hex", "[MTK],[Area_Integration],[Area_Hex8]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0, 0.0},
                                          {1.0, 0.0, 0.5},
                                          {1.0, 1.0, 0.5},
                                          {0.0, 1.0, 0.0},
                                          {0.0, 0.3, 1.0},
                                          {1.0, 0.3, 1.5},
                                          {1.0, 1.3, 1.5},
                                          {0.0, 1.3, 1.0}};
                Matrix<IdMat> tNodeLocalToGlobal = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify the local to global map
                //Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX8;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                std::string tFileOutput = "./mtk_hex8_cell_area_ut.exo";
                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // Setup Node Vertices (note: this data structure will be in the STK_Implementation
                moris::Cell<Vertex *> tElementVertices;
                for (size_t i = 0; i < tNodeIndices.numel(); i++)
                {
                    tElementVertices.push_back(new Vertex_STK(tNodeIds(i), tNodeIndices(i), tMesh));
                }

                // Setup cell associated with element index 0
                Cell_Info_Hex8 tHex8;
                Cell_STK tCell(&tHex8, 1, 0, tElementVertices, tMesh);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 8);
                CHECK(tCell.compute_cell_measure_general() == Approx(1.0));
            }
        }

    } // namespace mtk
} // namespace moris