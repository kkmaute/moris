#include "catch.hpp"
#include "math.h"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "fn_GEN_create_geometries.hpp"

#include "fn_GEN_check_equal.hpp"
#include "fn_GEN_create_simple_mesh.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> get_quad_local_coordinates(uint tNodeNumber)
        {
            real tRadians = M_PI * (5 + 2 * tNodeNumber) / 4;
            return sqrt(2) * Matrix<DDRMat>( {{cos(tRadians), sin(tRadians)}} );
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Linear Intersections", "[gen], [pdv], [intersection], [linear intersection]")
        {
            if (par_size() == 1)
            {
                // Create mesh
                mtk::Interpolation_Mesh* tMesh = nullptr;
                mtk::Integration_Mesh* tIGMesh = nullptr;
                create_simple_mesh(tMesh,tIGMesh, 2, 2);

                // Set up geometry
                Matrix<DDRMat> tADVs(0, 0);

                // Circle
                ParameterList tCircleParameterList = prm::create_geometry_parameter_list();
                tCircleParameterList.set("type", "circle");
                tCircleParameterList.set("constant_parameters", "-0.25, 0.0, 0.7499999999");
                tCircleParameterList.set("bspline_mesh_index", 0);

                // Plane
                ParameterList tPlaneParameterList = prm::create_geometry_parameter_list();
                tPlaneParameterList.set("type", "plane");
                tPlaneParameterList.set("constant_parameters", "0.25, 0.0, 1.0, 0.0");

                std::shared_ptr<mtk::Mesh_Manager> tMeshManager =
                        std::make_shared< mtk::Mesh_Manager >();

                tMeshManager->register_mesh_pair(tMesh, tIGMesh );

                // Create geometry engine
                Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries =
                        create_geometries({tCircleParameterList, tPlaneParameterList}, tADVs);
                tGeometryEngineParameters.mGeometries(0)->set_mesh( tMeshManager );
                tGeometryEngineParameters.mGeometries(0)->set_mesh_index( 0 );

                Geometry_Engine tGeometryEngine(tMesh, tGeometryEngineParameters);

                // TODO ensure this writes the mesh/fields correctly instead of just relying on no errors being thrown
                tGeometryEngine.output_fields_on_mesh(tMesh, "intersection_test.exo");

                // Solution for is_intersected() per geometry and per element
                Cell<Cell<bool>> tIsElementIntersected = {{true, true, true, true}, {false, true, false, true}};

                // Per geometry, per element, per edge
                Cell<Cell<Cell<bool>>> tIsEdgeIntersected = {{
                        {false, true, true, false},   // Geometry 0, Element 0
                        {false, false, true, true},   // Geometry 0, Element 1
                        {true, true, false, false},   // Geometry 0, Element 2
                        {true, false, false, true}},{ // Geometry 0, Element 3
                        {false, false, false, false}, // Geometry 1, Element 0
                        {true, false, true, false},   // Geometry 1, Element 1
                        {false, false, false, false}, // Geometry 1, Element 2
                        {true, false, true, false}}}; // Geometry 1, Element 3

                // Intersection coordinates
                real tFrac = 2.0 / (3.0 + sqrt(17.0));
                Matrix<DDRMat> tIntersectionLocalCoordinates = {{
                        -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5}};
                Cell<Matrix<DDRMat>> tIntersectionGlobalCoordinates = {
                        {{0.0, -0.5 - (tFrac / 2)}},
                        {{-1.0, 0.0}},
                        {{0.5, 0.0}},
                        {{0.0, -0.5 - (tFrac / 2)}},
                        {{-1.0, 0.0}},
                        {{0.0, 0.5 + (tFrac / 2)}},
                        {{0.5, 0.0}},
                        {{0.0, 0.5 + (tFrac / 2)}},
                        {{0.25, -1.0}},
                        {{0.25, 0.0}},
                        {{0.25, 0.0}},
                        {{0.25, 1.0}}};

                // Check element intersections
                uint tIntersectionCount = 0;
                for (uint tGeometryIndex = 0; tGeometryIndex < 2; tGeometryIndex++)
                {
                    for (uint tElementIndex = 0; tElementIndex < 4; tElementIndex++)
                    {
                        // Node indices per element
                        Matrix<IndexMat> tNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds(tElementIndex);

                        // Check edges for properly queued intersections
                        Matrix<DDRMat> tNodeCoordinates(4, 2);
                        for (uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++)
                        {
                            // Queue intersection
                            bool tIntersectionQueued = tGeometryEngine.queue_intersection( // TODO
                                    tNodeIndices(tNodeNumber),
                                    tNodeIndices((tNodeNumber + 1) % 4),
                                    {{}},
                                    {{}},
                                    tMesh->get_node_coordinate(tNodeIndices(tNodeNumber)),
                                    tMesh->get_node_coordinate(tNodeIndices((tNodeNumber + 1) % 4)),
                                    {{}},
                                    {});
                            REQUIRE(tIntersectionQueued == tIsEdgeIntersected(tGeometryIndex)(tElementIndex)(tNodeNumber));

                            // Check queued intersection
                            if (tIntersectionQueued)
                            {
                                // Check parents
                                bool tFirstParentOnInterface = false;
                                bool tSecondParentOnInterface = false;

                                // TODO when XTK supports background nodes on the interface, uncomment this
//                                if (tNodeIndices(tNodeNumber) == 3)
//                                {
//                                    tFirstParentOnInterface = true;
//                                }
//                                if (tNodeIndices((tNodeNumber + 1) % 4) == 3)
//                                {
//                                    tSecondParentOnInterface = true;
//                                }
                                CHECK(tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface);
                                CHECK(tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface);

                                // Check local coordinates
                                CHECK(tGeometryEngine.get_queued_intersection_local_coordinate() ==
                                        Approx(tIntersectionLocalCoordinates(tIntersectionCount)));

                                // Check global coordinates
                                check_equal(
                                        tGeometryEngine.get_queued_intersection_global_coordinates(),
                                        tIntersectionGlobalCoordinates(tIntersectionCount));

                                // Admit intersection
                                tGeometryEngine.admit_queued_intersection(9 + tIntersectionCount);

                                // Increment intersection count
                                tIntersectionCount++;
                            }

                            // Set node coordinates for element checking
                            tNodeCoordinates.set_row(tNodeNumber, tMesh->get_node_coordinate(tNodeIndices(tNodeNumber)));
                        }

                        // Check with solution
                        CHECK(tGeometryEngine.is_intersected(tNodeIndices, tNodeCoordinates) ==
                                tIsElementIntersected(tGeometryIndex)(tElementIndex));
                    }

                    // FIXME
                    if (tGeometryIndex == 0)
                    {
                        tGeometryEngine.advance_geometry_index();
                    }
                }

                // Check total number of intersections
                CHECK(tIntersectionCount == 12);

                // Test the new child nodes on the level set field (geometry 0)
                CHECK(tGeometryEngine.get_field_value(0,  9, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 10, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 11, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 12, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 13, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 14, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 15, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 16, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 17, {{}}) == Approx(0.423278));
                CHECK(tGeometryEngine.get_field_value(0, 18, {{}}) == Approx(-0.25));
                CHECK(tGeometryEngine.get_field_value(0, 19, {{}}) == Approx(-0.25));
                CHECK(tGeometryEngine.get_field_value(0, 20, {{}}) == Approx(0.423278));

                // Get the PDV host manager and set the number of total nodes
                Pdv_Host_Manager* tPdvHostManager = dynamic_cast<Pdv_Host_Manager*>(tGeometryEngine.get_design_variable_interface());

                // Test that the new intersections have been added to the PDV host manager, but ONLY for the circle
                Cell<Matrix<DDRMat>> tPdvValues(0);
                Cell<Matrix<DDSMat>> tIsActive(0);
                tPdvHostManager->get_ig_pdv_value(
                        {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}},
                        {PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE},
                        tPdvValues,
                        tIsActive);

                // Background nodes
                for (uint tNodeIndex = 0; tNodeIndex <= 8; tNodeIndex++)
                {
                    // Check if not active
                    CHECK(tIsActive(0)(tNodeIndex) == false);
                    CHECK(tIsActive(1)(tNodeIndex) == false);
                }

                // Nodes on the circle interface (depends on ADVs, active)
                for (uint tNodeIndex = 9; tNodeIndex <= 16; tNodeIndex++)
                {
                    // Check if active
                    CHECK(tIsActive(0)(tNodeIndex) == true);
                    CHECK(tIsActive(1)(tNodeIndex) == true);

                    // Check for node coordinates as PDV values
                    CHECK(tPdvValues(0)(tNodeIndex) == Approx(tIntersectionGlobalCoordinates(tNodeIndex - 9)(0)));
                    CHECK(tPdvValues(1)(tNodeIndex) == Approx(tIntersectionGlobalCoordinates(tNodeIndex - 9)(1)));
                }

                // Nodes on the plane interface (inactive) TODO Subset are still on circle
                for (uint tNodeIndex = 17; tNodeIndex <= 20; tNodeIndex++)
                {
                    // Check if not active
                    CHECK(tIsActive(0)(tNodeIndex) == false);
                    CHECK(tIsActive(1)(tNodeIndex) == false);
                }

                // TODO check sensitivities

                //------------------------------------------------------------------------------------------------------
                // Start second check
                //------------------------------------------------------------------------------------------------------

                // Set new ADVs, level set field now has no intersections
                tGeometryEngine.set_advs(Matrix<DDRMat>(16, 1, 1.0));
                tGeometryEngine.distribute_advs(tMesh);

                // Solution for is_intersected() per geometry and per element
                tIsElementIntersected = {{false, false, false, false}, {false, true, false, true}};

                // Per geometry, per element, per edge
                tIsEdgeIntersected = {{
                        {false, false, false, false},   // Geometry 0, Element 0
                        {false, false, false, false},   // Geometry 0, Element 1
                        {false, false, false, false},   // Geometry 0, Element 2
                        {false, false, false, false}},{ // Geometry 0, Element 3
                        {false, false, false, false},   // Geometry 1, Element 0
                        {true, false, true, false},     // Geometry 1, Element 1
                        {false, false, false, false},   // Geometry 1, Element 2
                        {true, false, true, false}}};   // Geometry 1, Element 3

                // Intersection local coordinates
                tIntersectionLocalCoordinates = {{-0.5, 0.5, -0.5, 0.5}};

                // Check element intersections
                tIntersectionCount = 0;
                for (uint tGeometryIndex = 0; tGeometryIndex < 2; tGeometryIndex++)
                {
                    for (uint tElementIndex = 0; tElementIndex < 4; tElementIndex++)
                    {
                        // Node indices per element
                        Matrix<IndexMat> tNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds(tElementIndex);

                        // Check edges for properly queued intersections
                        Matrix<DDRMat> tNodeCoordinates(4, 2);
                        for (uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++)
                        {
                            // Queue intersection
                            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                    tNodeIndices(tNodeNumber),
                                    tNodeIndices((tNodeNumber + 1) % 4),
                                    {{}},
                                    {{}},
                                    tMesh->get_node_coordinate(tNodeIndices(tNodeNumber)),
                                    tMesh->get_node_coordinate(tNodeIndices((tNodeNumber + 1) % 4)),
                                    {{}},
                                    {});
                            REQUIRE(tIntersectionQueued == tIsEdgeIntersected(tGeometryIndex)(tElementIndex)(tNodeNumber));

                            // Check queued intersection
                            if (tIntersectionQueued)
                            {
                                // Check parents
                                bool tFirstParentOnInterface = false;
                                bool tSecondParentOnInterface = false;

                                // TODO when XTK supports background nodes on the interface, uncomment this
//                                if (tNodeIndices(tNodeNumber) == 3)
//                                {
//                                    tFirstParentOnInterface = true;
//                                }
//                                if (tNodeIndices((tNodeNumber + 1) % 4) == 3)
//                                {
//                                    tSecondParentOnInterface = true;
//                                }
                                CHECK(tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface);
                                CHECK(tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface);

                                // Check local coordinates
                                CHECK(tGeometryEngine.get_queued_intersection_local_coordinate() ==
                                        Approx(tIntersectionLocalCoordinates(tIntersectionCount)));

                                // Admit intersection
                                tGeometryEngine.admit_queued_intersection(9 + tIntersectionCount);

                                // Increment intersection count
                                tIntersectionCount++;
                            }

                            // Set node coordinates for element checking
                            tNodeCoordinates.set_row(tNodeNumber, tMesh->get_node_coordinate(tNodeIndices(tNodeNumber)));
                        }

                        // Check with solution
                        CHECK(tGeometryEngine.is_intersected(tNodeIndices, tNodeCoordinates) ==
                                tIsElementIntersected(tGeometryIndex)(tElementIndex));
                    }

                    // FIXME
                    if (tGeometryIndex == 0)
                    {
                        tGeometryEngine.advance_geometry_index();
                    }
                }

                // Check total number of intersections
                CHECK(tIntersectionCount == 4);

                // Test the new child nodes on the level set field (geometry 0)
                CHECK(tGeometryEngine.get_field_value(0, 9,  {{}}) == Approx(1.0));
                CHECK(tGeometryEngine.get_field_value(0, 10, {{}}) == Approx(1.0));
                CHECK(tGeometryEngine.get_field_value(0, 11, {{}}) == Approx(1.0));
                CHECK(tGeometryEngine.get_field_value(0, 12, {{}}) == Approx(1.0));

                // Clean up
                delete tMesh;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Bilinear Intersections", "[gen], [pdv], [intersection], [bilinear intersection]")
        {
            if (par_size() == 1)
            {
                // Create mesh
                mtk::Interpolation_Mesh* tMesh = nullptr;
                mtk::Integration_Mesh* tIGMesh = nullptr;
                create_simple_mesh(tMesh,tIGMesh,2, 2);

                // Set up circle
                ParameterList tCircleParameterList = prm::create_geometry_parameter_list();
                tCircleParameterList.set("type", "circle");
                tCircleParameterList.set("constant_parameters", "-0.25, 0.0, 0.7499999999");
                tCircleParameterList.set("bspline_mesh_index", 0);
                tCircleParameterList.set("multilinear_intersections", true);
                Matrix<DDRMat> tADVs(0, 0);

                std::shared_ptr<mtk::Mesh_Manager> tMeshManager =
                          std::make_shared< mtk::Mesh_Manager >();

                  tMeshManager->register_mesh_pair(tMesh, tIGMesh );

                // Create geometry engine
                Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries = create_geometries({tCircleParameterList}, tADVs);
                tGeometryEngineParameters.mGeometries(0)->set_mesh( tMeshManager );
                tGeometryEngineParameters.mGeometries(0)->set_mesh_index( 0 );

                Geometry_Engine tGeometryEngine(tMesh, tGeometryEngineParameters);

                // Solution for is_intersected() per geometry and per element
                Cell<bool> tIsElementIntersected = {true, true, true, true};

                // Per element, per edge
                Cell<Cell<bool>> tIsEdgeIntersected = {
                        {false, true, true, false},  // Element 0
                        {false, false, true, true},  // Element 1
                        {true, true, false, false},  // Element 2
                        {true, false, false, true}}; // Element 3

                // Intersection coordinates
                real tFrac = 2.0 / (3.0 + sqrt(17.0));
                Matrix<DDRMat> tIntersectionLocalCoordinates = {{
                        -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5}};
                Cell<Matrix<DDRMat>> tIntersectionGlobalCoordinates = {
                        {{0.0, -0.5 - (tFrac / 2)}},
                        {{-1.0, 0.0}},
                        {{0.5, 0.0}},
                        {{0.0, -0.5 - (tFrac / 2)}},
                        {{-1.0, 0.0}},
                        {{0.0, 0.5 + (tFrac / 2)}},
                        {{0.5, 0.0}},
                        {{0.0, 0.5 + (tFrac / 2)}}};

                // Check element intersections
                uint tIntersectionCount = 0;
                for (uint tElementIndex = 0; tElementIndex < 4; tElementIndex++)
                {
                    // Get element info
                    Matrix<IndexMat> tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds(tElementIndex);
                    Matrix<DDUMat> tNodeIndices(4, 1);
                    Cell<Matrix<DDRMat>> tNodeCoordinates(4);
                    for (uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++)
                    {
                        tNodeIndices(tNodeNumber) = tSignedNodeIndices(tNodeNumber);
                        tNodeCoordinates(tNodeNumber) = tMesh->get_node_coordinate(tNodeIndices(tNodeNumber));
                    }

                    // Check edges for properly queued intersections
                    for (uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++)
                    {
                        // Queue intersection
                        bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                tNodeIndices(tNodeNumber),
                                tNodeIndices((tNodeNumber + 1) % 4),
                                get_quad_local_coordinates(tNodeNumber),
                                get_quad_local_coordinates((tNodeNumber + 1) % 4),
                                tNodeCoordinates(tNodeNumber),
                                tNodeCoordinates((tNodeNumber + 1) % 4),
                                tNodeIndices,
                                tNodeCoordinates);
                        REQUIRE(tIntersectionQueued == tIsEdgeIntersected(tElementIndex)(tNodeNumber));

                        // Check queued intersection
                        if (tIntersectionQueued)
                        {
                            // Check parents
                            bool tFirstParentOnInterface = false;
                            bool tSecondParentOnInterface = false;

                            // TODO when XTK supports background nodes on the interface, uncomment this
//                                if (tNodeIndices(tNodeNumber) == 3)
//                                {
//                                    tFirstParentOnInterface = true;
//                                }
//                                if (tNodeIndices((tNodeNumber + 1) % 4) == 3)
//                                {
//                                    tSecondParentOnInterface = true;
//                                }
                            CHECK(tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface);
                            CHECK(tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface);

                            // Check local coordinates
                            CHECK(tGeometryEngine.get_queued_intersection_local_coordinate() ==
                                    Approx(tIntersectionLocalCoordinates(tIntersectionCount)));

                            // Check global coordinates
                            CHECK(tGeometryEngine.get_queued_intersection_global_coordinates()(0) ==
                                    Approx(tIntersectionGlobalCoordinates(tIntersectionCount)(0)));
                            CHECK(tGeometryEngine.get_queued_intersection_global_coordinates()(1) ==
                                    Approx(tIntersectionGlobalCoordinates(tIntersectionCount)(1)));

                            // Admit intersection
                            tGeometryEngine.admit_queued_intersection(9 + tIntersectionCount);

                            // Increment intersection count
                            tIntersectionCount++;
                        }
                    }
                }

                // Check total number of intersections
                CHECK(tIntersectionCount == 8);

                // Test the new child nodes on the level set field
                CHECK(tGeometryEngine.get_field_value(0,  9, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 10, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 11, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 12, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 13, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 14, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 15, {{}}) == Approx(0.0));
                CHECK(tGeometryEngine.get_field_value(0, 16, {{}}) == Approx(0.0));

                // Get full element info for element 0
                Matrix<IndexMat> tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds(0);
                Matrix<DDUMat> tNodeIndices(4, 1);
                Cell<Matrix<DDRMat>> tNodeCoordinates(4);
                for (uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++)
                {
                    tNodeIndices(tNodeNumber) = tSignedNodeIndices(tNodeNumber);
                    tNodeCoordinates(tNodeNumber) = tMesh->get_node_coordinate(tNodeIndices(tNodeNumber));
                }

                // Queue custom intersection 1 and check for bilinear intersection
                bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                        0, 0, {{-1.0, -1.0}}, {{1.0, 1.0}}, {{}}, {{}}, tNodeIndices, tNodeCoordinates);
                REQUIRE(tIntersectionQueued == true);
                CHECK(tGeometryEngine.get_queued_intersection_local_coordinate() == Approx(0.137725));
                check_equal(tGeometryEngine.get_queued_intersection_global_coordinates(), {{-0.43113736, -0.43113736}});

                // Queue custom intersection 2 and check for no bilinear intersection
                tIntersectionQueued = tGeometryEngine.queue_intersection(
                        0, 0, {{1.0, -1.0}}, {{-1.0, 1.0}}, {{}}, {{}}, tNodeIndices, tNodeCoordinates);
                REQUIRE(tIntersectionQueued == false);

                // Queue custom intersection 3 and check for bilinear intersection
                tIntersectionQueued = tGeometryEngine.queue_intersection(
                        0, 0, {{0.75, 0.0}}, {{-0.75, 0.0}}, {{}}, {{}}, tNodeIndices, tNodeCoordinates);
                REQUIRE(tIntersectionQueued == true);
                CHECK(tGeometryEngine.get_queued_intersection_local_coordinate() == Approx(-0.520518));
                check_equal(tGeometryEngine.get_queued_intersection_global_coordinates(), {{-0.304805898, -0.5}});

                // Clean up
                delete tMesh;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
