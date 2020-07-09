#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometry.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_Exec_load_user_library.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{

    //------------------------------------------------------------------------------------------------------------------

    // Dummy user-defined geometry
    real user_defined_geometry_field(const moris::Matrix<DDRMat>& aCoordinates,
                                     const moris::Cell<real*>&    aParameters);

    void user_defined_geometry_sensitivity(const Matrix<DDRMat>& aCoordinates,
                                           const Cell<real*>&    aParameters,
                                                 Matrix<DDRMat>& aSensitivities);

    // Approximate check for DDRMat vector
    void check_approx(Matrix<DDRMat> aMat1, Matrix<DDRMat> aMat2);

    //------------------------------------------------------------------------------------------------------------------

    namespace ge
    {
        TEST_CASE("Circle Test", "[GE], [GE_CIRCLE]")
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
            Matrix<DDRMat> tADVs(5, 1);
            tADVs(0) = 0.0;
            tADVs(1) = 1.0;
            tADVs(2) = 2.0;
            tADVs(3) = 1.0;
            tADVs(4) = 2.0;
            std::shared_ptr<Geometry> tCircle1 = create_geometry(tCircle1ParameterList, tADVs);
            std::shared_ptr<Geometry> tCircle2 = create_geometry(tCircle2ParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0(1, 2, 0.0);
            Matrix<DDRMat> tCoordinates1(1, 2, 1.0);
            Matrix<DDRMat> tCoordinates2(1, 2, 2.0);

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

        TEST_CASE("Sphere Test", "[GE], [GE_SPHERE]")
        {
            // Set up geometry
            ParameterList tSphereParameterList = prm::create_geometry_parameter_list();
            tSphereParameterList.set("type", "sphere");
            tSphereParameterList.set("geometry_variable_indices", "all");
            tSphereParameterList.set("adv_indices", "all");

            // Create sphere
            Matrix<DDRMat> tADVs(4, 1);
            tADVs(0) = -1.0;
            tADVs(1) = 0.0;
            tADVs(2) = 1.0;
            tADVs(3) = 2.0;
            std::shared_ptr<Geometry> tSphere = create_geometry(tSphereParameterList, tADVs);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates0(1, 3, 0.0);
            Matrix<DDRMat> tCoordinates1(1, 3, 1.0);
            Matrix<DDRMat> tCoordinates2(1, 3, 2.0);

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
            tADVs(0) = 0.0;
            tADVs(3) = 1.0;
            tCoordinates1(2) = -1.0;
            tCoordinates2(1) = -2.0;

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

        TEST_CASE("User-Defined Geometry Test", "[GE], [GE_USER_DEFINED_GEOMETRY]")
        {
            // Create user-defined geometry
            Matrix<DDRMat> tADVs(2, 1);
            tADVs(0) = -1.0;
            tADVs(1) = 0.5;

            std::shared_ptr<Geometry> tGeometry = std::make_shared<User_Defined_Geometry>(tADVs,
                                                                                          Matrix<DDUMat>({{1, 0}}),
                                                                                          Matrix<DDUMat>({{0, 1}}),
                                                                                          Matrix<DDRMat>({{}}),
                                                                                          &user_defined_geometry_field,
                                                                                          &user_defined_geometry_sensitivity);

            // Set coordinates for checking
            Matrix<DDRMat> tCoordinates1(1, 2, 1.0);
            Matrix<DDRMat> tCoordinates2(1, 2, 2.0);

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
            tADVs(0) = 2.0;
            tCoordinates1(0) = 0.0;
            tCoordinates2(1) = -1.0;

            // Check field values
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates1) == Approx(8.0));
            CHECK(tGeometry->evaluate_field_value(0, tCoordinates2) == Approx(-7.5));

            // Check sensitivity values
            tGeometry->evaluate_sensitivity(0, tCoordinates1, tSensitivities);
            check_approx(tSensitivities, {{0.0, 12.0}});
            tGeometry->evaluate_sensitivity(0, tCoordinates2, tSensitivities);
            check_approx(tSensitivities, {{2.0, -12.0}});
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    real user_defined_geometry_field(const moris::Matrix<DDRMat>& aCoordinates,
                                     const moris::Cell<real*>&    aParameters)
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
