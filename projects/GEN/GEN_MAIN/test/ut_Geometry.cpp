#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_Exec_load_user_library.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_PRM_HMR_Parameters.hpp"

namespace moris
{

    //------------------------------------------------------------------------------------------------------------------

    // Dummy user-defined geometry
    real user_defined_geometry_field(const Matrix<DDRMat>& aCoordinates,
                                     const Cell<real*>&    aParameters);

    void user_defined_geometry_sensitivity(const Matrix<DDRMat>& aCoordinates,
                                           const Cell<real*>&    aParameters,
                                                 Matrix<DDRMat>& aSensitivities);

    // Approximate check for DDRMat vector
    void check_approx(Matrix<DDRMat> aMat1, Matrix<DDRMat> aMat2);

    //------------------------------------------------------------------------------------------------------------------

    namespace ge
    {

        // Check for ellipse location in a swiss cheese
        void check_swiss_cheese(std::shared_ptr<Geometry> aSwissCheese,
                                real aXCenter,
                                real aYCenter,
                                real aXSemidiameter,
                                real aYSemidiameter,
                                bool aCheck = true);

        TEST_CASE("Circle test", "[GEN], [geometry], [circle]")
        {
            // Set up geometry
            ParameterList tCircle1ParameterList = prm::create_geometry_parameter_list();
            tCircle1ParameterList.set("type", "circle");
            tCircle1ParameterList.set("geometry_variable_indices", "all");
            tCircle1ParameterList.set("adv_indices", "0, 1, 3");

            ParameterList tCircle2ParameterList = prm::create_geometry_parameter_list();
            tCircle2ParameterList.set("type", "circle");
            tCircle2ParameterList.set("geometry_variable_indices", "all");
            tCircle2ParameterList.set("adv_indices", "0, 2, 4");

            // Create circles
            Matrix<DDRMat> tADVs = {{0.0, 1.0, 2.0, 1.0, 2.0}};
            std::shared_ptr<Geometry> tCircle1 = create_geometry(tCircle1ParameterList, tADVs);
            std::shared_ptr<Geometry> tCircle2 = create_geometry(tCircle2ParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0 = {{0.0, 0.0}};
            Matrix<DDRMat> tCoordinates1 = {{1.0, 1.0}};
            Matrix<DDRMat> tCoordinates2 = {{2.0, 2.0}};

            // Check field values
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates1) == Approx(0.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(2.0) - 2.0));
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates2) == Approx(sqrt(5.0) - 1.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tCircle1->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 1.0, 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.0, 1.0, 0.0, -1.0}});
            tCircle1->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, -1.0}});
            tCircle1->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-2.0 / sqrt(5.0), -1.0 / sqrt(5.0), 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, 0.0, -1.0}});

            // Change ADVs and coordinates
            tADVs(0) = 1.0;
            tADVs(3) = 2.0;
            tADVs(4) = 3.0;
            tCoordinates0(0) = 1.0;
            tCoordinates0(1) = -1.0;
            tCoordinates1(0) = 3.0;
            tCoordinates1(1) = 1.0;
            tCoordinates2(0) = 4.0;
            tCoordinates2(1) = 2.0;

            // Check field values
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates1) == Approx(0.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(5.0) - 3.0));
            CHECK(tCircle1->evaluate_field_value(0, tCoordinates2) == Approx(sqrt(10.0) - 2.0));
            CHECK(tCircle2->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            tCircle1->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 1.0, 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.0, 1.0, 0.0, -1.0}});
            tCircle1->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-2.0 / sqrt(5.0), 0.0, 1.0 / sqrt(5.0), 0.0, -1.0}});
            tCircle1->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-3.0 / sqrt(10.0), -1.0 / sqrt(10.0), 0.0, -1.0, 0.0}});
            tCircle2->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, 0.0, -1.0}});
        }

        TEST_CASE("Superellipse test", "[GEN], [geometry], [superellipse]")
        {
            // Set up geometry
            ParameterList tSuperellipseParameterList = prm::create_geometry_parameter_list();
            tSuperellipseParameterList.set("type", "superellipse");
            tSuperellipseParameterList.set("geometry_variable_indices", "all");
            tSuperellipseParameterList.set("adv_indices", "all");

            // Create circles
            Matrix<DDRMat> tADVs = {{3.0, 4.0, 1.0, 2.0, 3.0}};
            std::shared_ptr<Geometry> tSuperellipse = create_geometry(tSuperellipseParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0 = {{2.0, 2.0}};
            Matrix<DDRMat> tCoordinates1 = {{3.0, 3.0}};
            Matrix<DDRMat> tCoordinates2 = {{4.0, 4.0}};

            // Check field values
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates0) == Approx(pow(2.0, 1.0/3.0) - 1.0));
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates1) == Approx(-0.5));
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tSuperellipse->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{pow(2.0, -2.0 / 3.0), pow(2.0, -5.0 / 3.0),
                    -pow(2.0, -2.0 / 3.0), -pow(2.0, -5.0 / 3.0), -0.0970335}});
            tSuperellipse->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.5, 0.0, -0.25, 0.0}});
            tSuperellipse->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, -1.0, 0.0, 0.0}});

            // Change ADVs and coordinates
            tADVs = {{2.0, 1.0, 4.0, 3.0, 4.0}};
            tCoordinates0 = {{-2.0, 1.0}};
            tCoordinates1 = {{0.0, 2.5}};
            tCoordinates2 = {{2.0, 5.0}};

            // Check field values
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates1) == Approx(pow(2.0, -0.75) - 1.0));
            CHECK(tSuperellipse->evaluate_field_value(0, tCoordinates2) == Approx(1.0 / 3.0));

            // Check sensitivity values
            tSuperellipse->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.25, 0.0, -0.25, 0.0, 0.0}});
            tSuperellipse->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{pow(2.0, 0.25) / 8.0, -pow(2.0, -0.75) / 3.0,
                    -pow(2.0, -0.75) / 8.0, -pow(2.0, -0.75) / 6.0, -0.0257572}});
            tSuperellipse->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{0.0, -1.0 / 3.0, 0.0, -4.0 / 9.0, 0.0}});
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Sphere test", "[GEN], [geometry], [sphere]")
        {
            // Set up geometry
            ParameterList tSphereParameterList = prm::create_geometry_parameter_list();
            tSphereParameterList.set("type", "sphere");
            tSphereParameterList.set("geometry_variable_indices", "all");
            tSphereParameterList.set("adv_indices", "all");

            // Create sphere
            Matrix<DDRMat> tADVs = {{-1.0, 0.0, 1.0, 2.0}};
            std::shared_ptr<Geometry> tSphere = create_geometry(tSphereParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0 = {{0.0, 0.0, 0.0}};
            Matrix<DDRMat> tCoordinates1 = {{1.0, 1.0, 1.0}};
            Matrix<DDRMat> tCoordinates2 = {{2.0, 2.0, 2.0}};

            // Check field values
            CHECK(tSphere->evaluate_field_value(0, tCoordinates0) == Approx(sqrt(2.0) - 2.0));
            CHECK(tSphere->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(5.0) - 2.0));
            CHECK(tSphere->evaluate_field_value(0, tCoordinates2) == Approx(sqrt(14.0) - 2.0));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tSphere->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{-sqrt(2.0)/ 2.0, 0.0, sqrt(2.0) / 2.0, -1.0}});
            tSphere->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-2.0 / sqrt(5.0), -1.0 / sqrt(5.0), 0.0, -1.0}});
            tSphere->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-3.0 / sqrt(14.0), -sqrt(2.0 / 7.0), -1.0 / sqrt(14.0), -1.0}});

            // Change ADVs and coordinates
            tADVs = {{0.0, 0.0, 1.0, 1.0}};
            tCoordinates1 = {{1.0, 1.0, -1.0}};
            tCoordinates2 = {{2.0, -2.0, 2.0}};

            // Check field values
            CHECK(tSphere->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tSphere->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(6.0) - 1.0));
            CHECK(tSphere->evaluate_field_value(0, tCoordinates2) == Approx(2.0));

            // Check sensitivity values
            tSphere->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.0, 1.0, -1.0}});
            tSphere->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-1.0 / sqrt(6.0), -1.0 / sqrt(6.0), sqrt(2.0 / 3.0), -1.0}});
            tSphere->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, -1.0}});
        }

        TEST_CASE("Superellipsoid test", "[GEN], [geometry], [superellipsoid]")
        {
            // Set up geometry
            ParameterList tSuperellipsoidParameterList = prm::create_geometry_parameter_list();
            tSuperellipsoidParameterList.set("type", "superellipsoid");
            tSuperellipsoidParameterList.set("geometry_variable_indices", "all");
            tSuperellipsoidParameterList.set("adv_indices", "all");

            // Create circles
            Matrix<DDRMat> tADVs = {{3.0, 4.0, 5.0, 1.0, 2.0, 4.0, 3.0}};
            std::shared_ptr<Geometry> tSuperellipsoid = create_geometry(tSuperellipsoidParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0 = {{2.0, 2.0, 5.0}};
            Matrix<DDRMat> tCoordinates1 = {{3.0, 3.0, 5.0}};
            Matrix<DDRMat> tCoordinates2 = {{4.0, 4.0, 5.0}};

            // Check field values
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates0) == Approx(pow(2.0, 1.0/3.0) - 1.0));
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates1) == Approx(-0.5));
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{pow(2.0, -2.0 / 3.0), pow(2.0, -5.0 / 3.0), 0.0,
                    -pow(2.0, -2.0 / 3.0), -pow(2.0, -5.0 / 3.0), 0.0, -0.0970335}});
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.5, 0.0, 0.0, -0.25, 0.0, 0.0}});
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0}});

            // Change ADVs and coordinates
            tADVs = {{2.0, 1.0, 0.0, 5.0, 4.0, 3.0, 4.0}};
            tCoordinates0 = {{2.0, -3.0, 0.0}};
            tCoordinates1 = {{2.0, -1.0, 1.5}};
            tCoordinates2 = {{2.0, 1.0, 4.0}};

            // Check field values
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates1) == Approx(pow(2.0, -0.75) - 1.0));
            CHECK(tSuperellipsoid->evaluate_field_value(0, tCoordinates2) == Approx(1.0 / 3.0));

            // Check sensitivity values
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.25, 0.0, 0.0, -0.25, 0.0, 0.0}});
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{0.0, pow(2.0, 0.25) / 8.0, -pow(2.0, -0.75) / 3.0,
                    0.0, -pow(2.0, -0.75) / 8.0, -pow(2.0, -0.75) / 6.0, -0.0257572}});
            tSuperellipsoid->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{0.0, 0.0, -1.0 / 3.0, 0.0, 0.0, -4.0 / 9.0, 0.0}});
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Level set geometry test", "[GEN], [geometry], [level set geometry]")
        {
            if (par_size() == 1)
            {
                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", "2, 2");
                tParameters.set( "domain_dimensions", "2, 2");
                tParameters.set( "domain_offset", "-1.0, -1.0");
                tParameters.set( "domain_sidesets", "1,2,3,4");
                tParameters.set( "lagrange_output_meshes", "0");

                tParameters.set( "lagrange_orders", "1");
                tParameters.set( "lagrange_pattern", "0");
                tParameters.set( "bspline_orders", "2");
                tParameters.set( "bspline_pattern", "0");
                
                tParameters.set( "lagrange_to_bspline", "0");
                
                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 3 );
                tParameters.set( "staircase_buffer", 3 );

                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement( 0 );
                tHMR.finalize();

                // Get interpolation mesh
                hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh(0);

                // Set up geometry
                real tRadius = 0.5;
                ParameterList tCircleParameterList = prm::create_geometry_parameter_list();
                tCircleParameterList.set("type", "circle");
                tCircleParameterList.set("constant_parameters", "0.0, 0.0, " + std::to_string(tRadius));
                tCircleParameterList.set("bspline_mesh_index", 0);
                Cell<std::shared_ptr<Geometry>> tGeometry(1);
                Matrix<DDRMat> tADVs(0, 0);
                tGeometry(0) = create_geometry(tCircleParameterList, tADVs);

                // Create geometry engine
                Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
                Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable, tInterpolationMesh);

                // Check values
                for (uint tNodeIndex = 0; tNodeIndex < tInterpolationMesh->get_num_nodes(); tNodeIndex++)
                {
                    Matrix<DDRMat> tNodeCoord = tInterpolationMesh->get_node_coordinate(tNodeIndex);
                    CHECK( tGeometryEngine.get_geometry_field_value(tNodeIndex, {{}}, 0) ==
                        Approx(sqrt(pow(tNodeCoord(0), 2) + pow(tNodeCoord(1), 2)) - tRadius) );
                }

                // Clean up
                delete tInterpolationMesh;

            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("User-defined geometry test", "[GEN], [geometry], [user-defined geometry]")
        {
            // Create user-defined geometry
            Matrix<DDRMat> tADVs = {{-1.0, 0.5}};

            std::shared_ptr<Geometry> tGeometry = std::make_shared<User_Defined_Geometry>(tADVs,
                                                                                          Matrix<DDUMat>({{1, 0}}),
                                                                                          Matrix<DDUMat>({{0, 1}}),
                                                                                          Matrix<DDRMat>({{}}),
                                                                                          &user_defined_geometry_field,
                                                                                          &user_defined_geometry_sensitivity);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates1 = {{1.0, 1.0}};
            Matrix<DDRMat> tCoordinates2 = {{2.0, 2.0}};

            // Check field values
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates1) == Approx(-0.75));
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates2) == Approx(-1.5));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tGeometry->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{1.0, 3.0}});
            tGeometry->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{2.0, 6.0}});

            // Change ADVs and coordinates
            tADVs = {{2.0, 0.5}};
            tCoordinates1 = {{0.0, 1.0}};
            tCoordinates2 = {{2.0, -1.0}};

            // Check field values
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates1) == Approx(8.0));
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates2) == Approx(-7.5));

            // Check sensitivity values
            tGeometry->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{0.0, 12.0}});
            tGeometry->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{2.0, -12.0}});
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Multigeometry test", "[GEN], [geometry], [multigeometry]")
        {
            // Set up geometry
            Cell<ParameterList> tCircleParameterLists(2);
            tCircleParameterLists(0) = prm::create_geometry_parameter_list();
            tCircleParameterLists(0).set("type", "circle");
            tCircleParameterLists(0).set("geometry_variable_indices", "all");
            tCircleParameterLists(0).set("adv_indices", "0, 1, 3");
            tCircleParameterLists(0).set("multigeometry_id", "circles");

            tCircleParameterLists(1) = prm::create_geometry_parameter_list();
            tCircleParameterLists(1).set("type", "circle");
            tCircleParameterLists(1).set("geometry_variable_indices", "all");
            tCircleParameterLists(1).set("adv_indices", "0, 2, 4");
            tCircleParameterLists(1).set("multigeometry_id", "circles");

            // Create multigeometry
            Matrix<DDRMat> tADVs = {{0.0, 1.0, 2.0, 1.0, 2.0}};
            Cell<std::shared_ptr<Geometry>> tGeometries = create_geometries(tCircleParameterLists, tADVs);

            // Should be only one total geometry
            REQUIRE(tGeometries.size() == 1);
            std::shared_ptr<Geometry> tMultigeometry = tGeometries(0);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0 = {{0.0, 0.0}};
            Matrix<DDRMat> tCoordinates1 = {{1.0, 1.0}};
            Matrix<DDRMat> tCoordinates2 = {{2.0, 2.0}};

            // Check field values
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(2.0) - 2.0));
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            Matrix<DDRMat> tSensitivities;
            tMultigeometry->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 1.0, 0.0, -1.0, 0.0}});
            tMultigeometry->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, -1.0}});
            tMultigeometry->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, 0.0, -1.0}});

            // Change ADVs and coordinates
            tADVs(0) = 1.0;
            tADVs(3) = 2.0;
            tADVs(4) = 3.0;
            tCoordinates0(0) = 1.0;
            tCoordinates0(1) = -1.0;
            tCoordinates1(0) = 3.0;
            tCoordinates1(1) = 1.0;
            tCoordinates2(0) = 4.0;
            tCoordinates2(1) = 2.0;

            // Check field values
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates0) == Approx(0.0));
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates1) == Approx(sqrt(5.0) - 3.0));
            CHECK(tMultigeometry->evaluate_field_value(0, tCoordinates2) == Approx(0.0));

            // Check sensitivity values
            tMultigeometry->evaluate_sensitivity(0, tCoordinates0, tSensitivities);
            check_approx(tSensitivities, {{0.0, 1.0, 0.0, -1.0, 0.0}});
            tMultigeometry->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{-2.0 / sqrt(5.0), 0.0, 1.0 / sqrt(5.0), 0.0, -1.0}});
            tMultigeometry->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{-1.0, 0.0, 0.0, 0.0, -1.0}});
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Swiss cheese slice test", "[GEN], [geometry], [swiss cheese slice]")
        {
            // Set up geometry
            ParameterList tSwissCheeseParameterList = prm::create_swiss_cheese_slice_parameter_list();
            tSwissCheeseParameterList.set("left_bound", -2.0);
            tSwissCheeseParameterList.set("right_bound", 2.0);
            tSwissCheeseParameterList.set("bottom_bound", -1.0);
            tSwissCheeseParameterList.set("top_bound", 1.0);
            tSwissCheeseParameterList.set("hole_x_semidiameter", 0.2);
            tSwissCheeseParameterList.set("hole_y_semidiameter", 0.1);
            tSwissCheeseParameterList.set("number_of_x_holes", 3);
            tSwissCheeseParameterList.set("number_of_y_holes", 5);

            // Create swiss cheese
            Matrix<DDRMat> tADVs = {{}};
            std::shared_ptr<Geometry> tSwissCheese = create_geometry(tSwissCheeseParameterList, tADVs);

            // Check holes
            check_swiss_cheese(tSwissCheese, -2.0, -1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, -2.0, 0.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, -2.0, 1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, -1.0, -1.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, -1.0, 0.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, -1.0, 1.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 0.0, -1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 0.0, 0.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 0.0, 1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 1.0, -1.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 1.0, 0.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 1.0, 1.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 2.0, -1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 2.0, 0.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 2.0, 1.0, 0.2, 0.1);
            check_swiss_cheese(tSwissCheese, 3.0, -1.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 3.0, 0.0, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.2, 0.1, false);
            check_swiss_cheese(tSwissCheese, 3.0, 1.0, 0.2, 0.1, false);

            // Create swiss cheese
            tSwissCheeseParameterList.set("number_of_x_holes", 0);
            tSwissCheeseParameterList.set("number_of_y_holes", 0);
            tSwissCheeseParameterList.set("hole_x_semidiameter", 0.1);
            tSwissCheeseParameterList.set("hole_y_semidiameter", 0.2);
            tSwissCheeseParameterList.set("target_x_spacing", 2.1);
            tSwissCheeseParameterList.set("target_y_spacing", 0.55);
            tSwissCheeseParameterList.set("row_offset", 1.0);
            tSwissCheese = create_geometry(tSwissCheeseParameterList, tADVs);
            
            // Check holes
            check_swiss_cheese(tSwissCheese, -2.0, -1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, -2.0, 0.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, -2.0, 1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, -1.0, -1.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, -1.0, 0.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, -1.0, 1.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 0.0, -1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 0.0, 0.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 0.0, 1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 1.0, -1.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 1.0, 0.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 1.0, 1.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 2.0, -1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 2.0, 0.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 2.0, 1.0, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 3.0, -1.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 3.0, 0.0, 0.1, 0.2, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.1, 0.2);
            check_swiss_cheese(tSwissCheese, 3.0, 1.0, 0.1, 0.2, false);
            
            // Create swiss cheese
            tSwissCheeseParameterList.set("hole_x_semidiameter", 0.45);
            tSwissCheeseParameterList.set("hole_y_semidiameter", 0.45);
            tSwissCheeseParameterList.set("target_y_spacing", 1.1);
            tSwissCheeseParameterList.set("allow_less_than_target_spacing", false);
            tSwissCheeseParameterList.set("row_offset", 0.0);
            tSwissCheese = create_geometry(tSwissCheeseParameterList, tADVs);
            
            // Check holes
            check_swiss_cheese(tSwissCheese, -2.0, -1.0, 0.45, 0.45);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -2.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -2.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -2.0, 1.0, 0.45, 0.45);
            check_swiss_cheese(tSwissCheese, -1.0, -1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -1.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -1.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, -1.0, 1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 0.0, -1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 0.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 0.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 0.0, 1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 1.0, -1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 1.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 1.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 1.0, 1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 2.0, -1.0, 0.45, 0.45);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 2.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 2.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 2.0, 1.0, 0.45, 0.45);
            check_swiss_cheese(tSwissCheese, 3.0, -1.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 3.0, 0.0, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 3.0, -0.5, 0.45, 0.45, false);
            check_swiss_cheese(tSwissCheese, 3.0, 1.0, 0.45, 0.45, false);
        }

        //------------------------------------------------------------------------------------------------------------------

        void check_swiss_cheese(std::shared_ptr<Geometry> aSwissCheese,
                                real aXCenter,
                                real aYCenter,
                                real aXSemidiameter,
                                real aYSemidiameter,
                                bool aCheck)
        {
            CHECK((aSwissCheese->evaluate_field_value(0, {{aXCenter + aXSemidiameter, aYCenter}}) == Approx(0.0)) == aCheck);
            CHECK((aSwissCheese->evaluate_field_value(0, {{aXCenter, aYCenter + aYSemidiameter}}) == Approx(0.0)) == aCheck);
            CHECK((aSwissCheese->evaluate_field_value(0, {{aXCenter - aXSemidiameter, aYCenter}}) == Approx(0.0)) == aCheck);
            CHECK((aSwissCheese->evaluate_field_value(0, {{aXCenter, aYCenter - aYSemidiameter}}) == Approx(0.0)) == aCheck);
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    real user_defined_geometry_field(const Matrix<DDRMat>& aCoordinates,
                                     const Cell<real*>&    aParameters)
    {
        return aCoordinates(0) * pow(*aParameters(0), 2) + aCoordinates(1) * pow(*aParameters(1), 3);
    }

    //------------------------------------------------------------------------------------------------------------------

    void user_defined_geometry_sensitivity(const Matrix<DDRMat>& aCoordinates,
                                           const Cell<real*>&    aParameters,
                                           Matrix<DDRMat>& aSensitivities)
    {
        aSensitivities = {{2 * aCoordinates(0) * *aParameters(0), 3 * aCoordinates(1) * pow(*aParameters(1), 2)}};
    }

    //------------------------------------------------------------------------------------------------------------------

    void check_approx(Matrix<DDRMat> aMat1, Matrix<DDRMat> aMat2)
    {
        REQUIRE(aMat1.length() == aMat2.length());
        for (uint tIndex = 0; tIndex < aMat1.length(); tIndex++)
        {
            CHECK(aMat1(tIndex) == Approx(aMat2(tIndex)));
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}
