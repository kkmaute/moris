//
// Created by christopherson on 5/19/20.
//

#include "catch.hpp"
#include "cl_Matrix.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_geometry.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_Exec_load_user_library.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"


namespace moris
{
    namespace ge
    {
        TEST_CASE("Discrete property based on ADVs", "[GE], [GE_DISCRETE_PROPERTY]")
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
            tParameterLists(2)(0).set("type", "discrete");
            tParameterLists(2)(0).set("name", "density");
            tParameterLists(2)(0).set("property_variable_indices", "all");
            tParameterLists(2)(0).set("adv_indices", "all");
            tParameterLists(2)(0).set("pdv_type", "DENSITY");

            Geometry_Engine tGeometryEngine(tParameterLists);
        }

        TEST_CASE("Property dependency test", "[GE], [GE_PROPERTY_DEPENDENCY]")
        {
            // Set up default parameter lists
            moris::Cell<moris::Cell<ParameterList>> tParameterLists(3);
            tParameterLists(0).resize(1);
            tParameterLists(2).resize(4);
            tParameterLists(0)(0) = moris::prm::create_gen_parameter_list();
            tParameterLists(2)(0) = moris::prm::create_gen_property_parameter_list();
            tParameterLists(2)(1) = moris::prm::create_gen_property_parameter_list();
            tParameterLists(2)(2) = moris::prm::create_gen_property_parameter_list();
            tParameterLists(2)(3) = moris::prm::create_gen_property_parameter_list();

            // Modify parameters
            tParameterLists(2)(0).set("type", "");
            tParameterLists(2)(0).set("name", "property_1");
            tParameterLists(2)(0).set("dependencies", "");

            tParameterLists(2)(1).set("type", "");
            tParameterLists(2)(1).set("name", "property_2");
            tParameterLists(2)(1).set("dependencies", "property_1");

            tParameterLists(2)(2).set("type", "");
            tParameterLists(2)(2).set("name", "property_3");
            tParameterLists(2)(2).set("dependencies", "property_1");

            tParameterLists(2)(3).set("type", "");
            tParameterLists(2)(3).set("name", "property_4");
            tParameterLists(2)(3).set("dependencies", "property_2,property_3");

            //Geometry_Engine tGeometryEngine(tParameterLists);
        }
    }   // end ge namespace
}   // end moris namespace
