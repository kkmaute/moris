#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometries.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_Exec_load_user_library.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Scaled_Field.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Discrete property based on ADVs", "[gen], [property], [discrete property]")
        {
            // Set up default parameter lists
            moris::Cell<moris::Cell<ParameterList>> tParameterLists(3);
            tParameterLists(0).resize(1);
            tParameterLists(2).resize(1);
            tParameterLists(0)(0) = moris::prm::create_gen_parameter_list();
            tParameterLists(2)(0) = moris::prm::create_gen_property_parameter_list();

            // Modify parameters
            tParameterLists(0)(0).set("advs_size", 4);
            tParameterLists(0)(0).set("initial_advs_fill", 0.3);
            tParameterLists(0)(0).set("lower_bounds_fill", 0.0);
            tParameterLists(0)(0).set("upper_bounds_fill", 1.0);
            tParameterLists(0)(0).set("PDV_types", "DENSITY");

            tParameterLists(2)(0).set("type", "discrete");
            tParameterLists(2)(0).set("name", "density");
            tParameterLists(2)(0).set("property_variable_indices", "all");
            tParameterLists(2)(0).set("adv_indices", "all");
            tParameterLists(2)(0).set("pdv_type", "DENSITY");

            Geometry_Engine tGeometryEngine(tParameterLists);
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE("Scaled field property", "[gen], [property], [scaled field]")
        {
            // Circle
            std::shared_ptr<Geometry> tCircle = std::make_shared<Circle>(0.0, 0.0, 0.5, "circle");

            // ADVs
            Matrix<DDRMat> tADVs(0, 0);

            // Random distribution
            std::uniform_real_distribution<real> tUniform(-100.0, 100.0);
            std::default_random_engine tEngine;

            for (uint tScaleRun = 0; tScaleRun < 4; tScaleRun++)
            {
                // Create scaled field
                real tScale = (real)tScaleRun;
                Scaled_Field tScaledField(tADVs, Matrix<DDUMat>(0, 0), Matrix<DDUMat>(0, 0), Matrix<DDRMat>(1, 1, tScale), {tCircle});

                for (uint tCoordinateCheck = 0; tCoordinateCheck < 4; tCoordinateCheck++)
                {
                    // Get random coordinates
                    Matrix<DDRMat> tCoordinates({{tUniform(tEngine), tUniform(tEngine)}});

                    // Checks
                    CHECK(tScaledField.get_field_value(0, tCoordinates) == Approx(tCircle->get_field_value(0, tCoordinates) * tScale));
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
