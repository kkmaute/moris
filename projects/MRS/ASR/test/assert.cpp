/*
 * assert.cpp
 *
 *  Created on: May 8, 2017
 *      Author: doble
 */

// third party
#include <catch.hpp>
#include "assert.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
namespace moris
{
TEST_CASE(
        "MORIS::assert",
        "[MORIS],[assert]")
                                {
    moris::uint p_size = 1;
    moris::uint p_rank = 1;
#ifdef MORIS_HAVE_PARALLEL
    p_size = par_size();
    p_rank = par_rank();
#endif

    SECTION("MORIS Assert Output Format Check")
    {
        // Turn this on if you want to check the format of the output message when a moris assert fails
        // No requirement tests (could add require throw)
        bool tOnOff = 0;
        if(tOnOff)
        {
            MORIS_ASSERT( false, "Test output message for MORIS ASSERT" );
        }
    }


    if(p_size>0)
        SECTION("Controlled MPI exit test")
        {
        // Processor 1 breaks the code and the other processors are expected to terminate the program when the receive the signal
        //
        bool tOnOff = 0;

#ifdef MORIS_USE_DEBUG
        if(p_rank == 1)
        {
            MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
        }
#endif

        if (tOnOff)
        {
            if(p_rank == 1)
            {
                MORIS_ASSERT( false, "Test output message for MORIS ASSERT" );
            }

            else
            {
                MORIS_ASSERT( true, "Test output message for MORIS ASSERT" );
            }
        }
        }
                                }


}

