#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometry.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        TEST_CASE("Geometry_Analytic with ADVs test", "[GE], [GE_GEOM_ADV]")
        {
            // Test parameter list
            // Set up default parameter lists
            moris::Cell<moris::Cell<ParameterList>> tParameterLists(2);
            tParameterLists(0).resize(1);
            tParameterLists(1).resize(2);
            tParameterLists(0)(0) = moris::prm::create_gen_parameter_list();
            tParameterLists(1)(0) = moris::prm::create_geometry_parameter_list();
            tParameterLists(1)(1) = moris::prm::create_geometry_parameter_list();

            // Modify parameters
            tParameterLists(0)(0).set("initial_advs", "0.0, 1.0, 2.0");

            tParameterLists(1)(0).set("type", "circle");
            tParameterLists(1)(0).set("geometry_variable_indices", "0, 2");
            tParameterLists(1)(0).set("adv_indices", "0, 1");
            tParameterLists(1)(0).set("constant_parameters", "1.0");

            tParameterLists(1)(1).set("type", "circle");
            tParameterLists(1)(1).set("geometry_variable_indices", "0, 2");
            tParameterLists(1)(1).set("adv_indices", "0, 2");
            tParameterLists(1)(1).set("constant_parameters", "2.0");

            // Create geometry separately
            Matrix<DDRMat> tADVs(3, 1);
            tADVs(0) = 0.0;
            tADVs(1) = 1.0;
            tADVs(2) = 2.0;
            std::shared_ptr<Geometry_Analytic> tCircle1 = create_geometry(tParameterLists(1)(0), tADVs);
            std::shared_ptr<Geometry_Analytic> tCircle2 = create_geometry(tParameterLists(1)(1), tADVs);

            // Check field values
            Matrix<DDRMat> tCoordinates(1, 2);
            tCoordinates(0) = 0.0;
            tCoordinates(1) = 0.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == 0.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == 0.0);
            tCoordinates(0) = 1.0;
            tCoordinates(1) = 1.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == 0.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == sqrt(2.0) - 2.0);
            tCoordinates(0) = 2.0;
            tCoordinates(1) = 2.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == sqrt(5.0) - 1.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == 0.0);

            // Change ADVs and check again
            tADVs(0) = 1.0;
            tADVs(1) = 2.0;
            tADVs(2) = 3.0;

            tCoordinates(0) = 1.0;
            tCoordinates(1) = -1.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == 0.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == 0.0);
            tCoordinates(0) = 3.0;
            tCoordinates(1) = 1.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == 0.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == sqrt(5.0) - 3.0);
            tCoordinates(0) = 4.0;
            tCoordinates(1) = 2.0;
            CHECK(tCircle1->evaluate_field_value(tCoordinates) == sqrt(10.0) - 2.0);
            CHECK(tCircle2->evaluate_field_value(tCoordinates) == 0.0);

        }   // end test case
    }   // end ge namespace
}   // end moris namespace
