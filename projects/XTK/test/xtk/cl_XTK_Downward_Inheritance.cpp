/*
 * cl_XTK_Downward_Inheritance.cpp
 *
 *  Created on: Jul 21, 2017
 *      Author: ktdoble
 */

#include <utility>

#include "containers/cl_XTK_Cell.hpp"
#include "xtk/cl_XTK_Downward_Inheritance.hpp"
#include "catch.hpp"

TEST_CASE("Downward Inheritance",
          "[XTK][INHERITANCE]")
        {
    /*
     * Tests the Downward inheritance structure to see if pairs are registered and stored correctly
     * and that existing pairs are not overwritten
     */
    xtk::Downward_Inheritance<xtk::size_t,xtk::size_t> tInheritance(5);
    tInheritance.register_new_inheritance_pair(1, 6);
    tInheritance.register_new_inheritance_pair(1, 7);
    CHECK(tInheritance.has_inheritance(1));
    CHECK(!tInheritance.has_inheritance(2));
    CHECK(tInheritance.get_inheritance(1)==6); // Make sure 7 did not overwrite 6

        }
