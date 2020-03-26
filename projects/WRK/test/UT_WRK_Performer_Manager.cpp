

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
	    std::string tStringMoris = std::getenv( "MORISROOT" );
	    std::string tString = tStringMoris + "/projects/FEM/MDL/test/data/Input_test.so";
    std::shared_ptr< Library_IO >tLibrary = std::make_shared< Library_IO >( tString );

	wrk::Performer_Manager tPerformerManager( tLibrary );

	tPerformerManager.initialize();

	tPerformerManager.perform();
    }
}
