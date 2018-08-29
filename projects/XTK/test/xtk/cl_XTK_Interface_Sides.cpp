/*
 * cl_XTK_Interface_Sides.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: ktdoble
 */




#include "catch.hpp"

#include<memory>

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg_typedefs.hpp"

// XTKL: Cell Container
#include "containers/cl_XTK_Cell.hpp"
#include "xtk/cl_XTK_Interface_Sides.hpp"

namespace xtk
{
TEST_CASE("Interface Side Container Test","[INTERFACE_SIDES]")
{

    Interface_Sides<size_t> tInterfaceSides(5);

    // create a side indices
    moris::Matrix<size_t,Default_Matrix_Integer> tSides({{4,5,7,1}});

    tInterfaceSides.add_interface_side_with_side_index(tSides);

    size_t tNumSides = tInterfaceSides.get_num_interface_sides();
    CHECK(tNumSides == 4);
    Cell<size_t> const & tSideCell = tInterfaceSides.get_interface_sides();

    CHECK(tSideCell(0) == 4);
    CHECK(tSideCell(1) == 5);
    CHECK(tSideCell(2) == 7);
    CHECK(tSideCell(3) == 1);

}
}
