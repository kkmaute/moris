#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "fn_GEN_create_properties.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Scaled_Field.hpp"
#include "fn_GEN_check_equal.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Constant property", "[gen], [property], [constant property]")
        {
            // Create ADVs
            Matrix<DDRMat> tADVs = {{0.0}};

            // Set up property
            ParameterList tConstantPropertyParameterList = prm::create_gen_property_parameter_list();
            tConstantPropertyParameterList.set("type", "constant");
            tConstantPropertyParameterList.set("field_variable_indices", "0");
            tConstantPropertyParameterList.set("adv_indices", "0");
            std::shared_ptr<Property> tConstantProperty = create_property(tConstantPropertyParameterList, tADVs);

            // Random distribution
            std::uniform_real_distribution<real> tUniform(-100.0, 100.0);
            std::default_random_engine tEngine;

            for (uint tTestRun = 0; tTestRun < 4; tTestRun++)
            {
                // Create scaled field
                tADVs(0) = tUniform(tEngine);

                // Loop over coordinate checks
                for (uint tCoordinateCheck = 0; tCoordinateCheck < 4; tCoordinateCheck++)
                {
                    // Get random coordinates
                    Matrix<DDRMat> tCoordinates({{tUniform(tEngine), tUniform(tEngine)}});

                    // Checks
                    CHECK(tConstantProperty->get_field_value(0, tCoordinates) == Approx(tADVs(0)));
                    check_equal(tConstantProperty->get_field_sensitivities(0, tCoordinates), {{1.0}});
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Scaled field property", "[gen], [property], [scaled field]")
        {
            // Create ADVs
            Matrix<DDRMat> tADVs = {{0.0, 0.0, 0.5}};

            // Set up and create geometry
            ParameterList tCircleParameterList = prm::create_geometry_parameter_list();
            tCircleParameterList.set("type", "circle");
            tCircleParameterList.set("name", "My Circle");
            tCircleParameterList.set("field_variable_indices", "0, 1, 2");
            tCircleParameterList.set("adv_indices", "0, 1, 2");
            Cell<std::shared_ptr<Geometry>> tGeometries = create_geometries({tCircleParameterList}, tADVs);
            std::shared_ptr<Geometry> tCircle = tGeometries(0);

            // Set up property
            ParameterList tScaledFieldParameterList = prm::create_gen_property_parameter_list();
            tScaledFieldParameterList.set("type", "scaled_field");
            tScaledFieldParameterList.set("dependencies", "My Circle");

            // Random distribution
            std::uniform_real_distribution<real> tUniform(-100.0, 100.0);
            std::default_random_engine tEngine;

            for (uint tTestRun = 0; tTestRun < 4; tTestRun++)
            {
                // Create scaled field
                real tScale = tUniform(tEngine);
                tScaledFieldParameterList.set("constant_parameters", std::to_string(tScale));
                Cell<std::shared_ptr<Property>> tProperties = create_properties({tScaledFieldParameterList}, tADVs, tGeometries);

                // Check that one property was created, and assign it as scaled field
                REQUIRE(tProperties.size() == 1);
                std::shared_ptr<Property> tScaledField = tProperties(0);

                // Loop over coordinate checks
                for (uint tCoordinateCheck = 0; tCoordinateCheck < 4; tCoordinateCheck++)
                {
                    // Get random coordinates
                    Matrix<DDRMat> tCoordinates({{tUniform(tEngine), tUniform(tEngine)}});

                    // Checks
                    CHECK(tScaledField->get_field_value(0, tCoordinates) ==
                            Approx(tCircle->get_field_value(0, tCoordinates) * tScale));
                    check_equal(tScaledField->get_field_sensitivities(0, tCoordinates),
                            Matrix<DDRMat>(tCircle->get_field_sensitivities(0, tCoordinates) * tScale));
                    check_equal(tScaledField->get_determining_adv_ids(0, tCoordinates),
                            tCircle->get_determining_adv_ids(0, tCoordinates));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
