/*
 * UT_Parameter_List.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: sonne
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

#include "paths.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_PRM_GEN_Parameters.hpp"
#include "fn_Exec_load_user_library.hpp"

//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
//------------------------------------------------------------------------------
TEST_CASE("property list test","[GE],[propListTest_00]")
{
    std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/libbenchmark04.so";

    std::shared_ptr< Library_IO > tLibrary = std::make_shared< Library_IO >( tInputFilePath );

    // load the MSI parameter list
    std::string tGENString = "GENParameterList";
    MORIS_PARAMETER_FUNCTION tGENParameterListFunc = tLibrary->load_parameter_file( tGENString );

    moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
    tGENParameterListFunc( tGENParameterList );

    GEN_Geometry_Engine tGE( tGENParameterList(0)(0) );

    tGE.initialize( tLibrary );

}   // end test case



}   // end ge namespace
}   // end moris namespace
