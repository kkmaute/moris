

#include "catch.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "fn_Exec_load_user_library.hpp"

using namespace moris;


TEST_CASE( "WRK_Test ", "[moris],[WRK_Test]" )
{
	if( par_size() == 1 )
	{
    std::shared_ptr< Library_IO >tLibrary = std::make_shared< Library_IO >( "/home/schmidt/codes/moris/projects/FEM/MDL/test/data/Input_test.so" );

	wrk::Performer_Manager tPerformerManager( tLibrary );

	tPerformerManager.initialize();

	tPerformerManager.perform();
    }
}
