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

//    // User-defined property functions
//    real property1_evaluate_field_value(const moris::Matrix< DDRMat >    & aCoordinates,
//                                     const moris::Cell< moris::real* > & aParameters)
//    {
//        return 1.0;
//    }
//    real property2_evaluate_field_value(const moris::Matrix< DDRMat >    & aCoordinates,
//                                        const moris::Cell< moris::real* > & aParameters)
//    {
//        return 1.0;
//    }
//    real property3_evaluate_field_value(const moris::Matrix< DDRMat >    & aCoordinates,
//                                        const moris::Cell< moris::real* > & aParameters)
//    {
//        return 1.0;
//    }
//    real property4_evaluate_field_value(const moris::Matrix< DDRMat >    & aCoordinates,
//                                        const moris::Cell< moris::real* > & aParameters)
//    {
//        return 1.0;
//    }
//
//    void property_evaluate_sensitivity(const moris::Matrix< DDRMat >    & aCoordinates,
//                                     const moris::Cell< moris::real* > & aParameters,
//                                     moris::Matrix< DDRMat >    & aSensitivities)
//    {
//        MORIS_ERROR(false, "sensitivities not defined");
//    }

    namespace ge
    {
        TEST_CASE("Property dependency test", "[GE], [GE_PROP]")
        {
            // Test parameter list
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

            tParameterLists(2)(0).set("type", "");
            tParameterLists(2)(0).set("name", "property_2");
            tParameterLists(2)(0).set("dependencies", "property_1");

            tParameterLists(2)(0).set("type", "");
            tParameterLists(2)(0).set("name", "property_3");
            tParameterLists(2)(0).set("dependencies", "property_1");

            tParameterLists(2)(0).set("type", "");
            tParameterLists(2)(0).set("name", "property_4");
            tParameterLists(2)(0).set("dependencies", "property_2,property_3");

            //Geometry_Engine tGeometryEngine(tParameterLists);
        }
    }   // end ge namespace
}   // end moris namespace
